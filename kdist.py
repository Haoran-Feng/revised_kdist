#
#     Calculates kinematic distances using the revised prescription given in
#     Reid et al 2009 (Paper VI):
#     Version 2 (projects Vlsr to Galactic Plane)
#     Version 3 (fixes bug: used Theta instead of Theta_0 in rootterm)
#     translated into IDL by Z. Jiang 2019.06.19
#     translated into Python by H. Feng 2019.07

import argparse
import pathlib
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import numpy as np
import json
import warnings
warnings.filterwarnings('ignore')


def parse_args():
    """
    Get command line arguments.
    :return:
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "source_file",
        type=str,
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help='Output filename.'
    )

    parser.add_argument(
        "--parameter_file",
        type=str,
        default="parameter_file.json"
    )

    parser.add_argument(
        "--use_radec",
        action="store_true"
    )

    parser.add_argument(
        "--max_iter",
        type=int,
        default=1000
    )

    parser.add_argument(
        '--comment_symbol',
        type=str,
        default='!'
    )

    return parser.parse_args()


def get_parameters_from_json(param_json_file: str) -> dict:
    """
    Get parameters from json file.
    :param param_json_file: str, path to the json file
    :return: dict of astropy Quantity, param names as keys, Quantities as values.
    """
    path = pathlib.Path(param_json_file)
    if not path.is_file():
        raise FileNotFoundError("Please check parameter file.")

    with open(param_json_file, 'r') as fp:
        params = json.load(fp)

    p_dict = dict()
    for k, v in params.items():
        value, unit, *_ = v
        p_dict[k] = value * u.Unit(unit)

    return p_dict


def get_source_data_from_file(src_file: str, if_use_ra_dec: bool, data_file_comment_symbol='!'):
    """
    Get source data from a text file
    :param src_file:  str, path to the file
    :param if_use_ra_dec: bool, set to True if ra dec coordinates are used in the data file
    :param data_file_comment_symbol: str, the symbol that denotes comment lines.
    :return: src_list (numpy.ndarray of source names),
             coords (astropy SkyCoord object),
             v_lsr (astropy Qantity array)
    """
    data = pd.read_csv(src_file, header=None, delim_whitespace=True, comment=data_file_comment_symbol)

    if if_use_ra_dec:
        ra = data.iloc[:, 1]
        dec = data.iloc[:, 2]
        coords = SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg')
        coords = coords.galactic

    else:
        l = data.iloc[:, 1]
        b = data.iloc[:, 2]
        coords = SkyCoord(l=l, b=b, frame='galactic', unit='deg')

    v_lsr = data.iloc[:, 3].values * u.Unit('km/s')
    farnear = data.iloc[:, 4].values
    src_list = data.iloc[:, 0].values

    # TODO: handle exceptions

    return src_list, coords, v_lsr, farnear


class DkCalculator(object):

    def __init__(self, params: dict):
        """
        :param params: dict that contains all the params. It MUST contains key-value pairs like this:
              {
                'Ro':         <Quantity 8.5 kpc>,
                'To':         <Quantity 220. km / s>,
                'dT/dr':      <Quantity -0.2 km / (kpc s)>,
                'Uo':         <Quantity 10.7 km / s>,
                'Vo':         <Quantity 15.6 km / s>,
                'Wo':         <Quantity 8.9 km / s>}
                'Ubar':       <Quantity 2.9 km / s>,
                'Vbar':       <Quantity -1.6 km / s>,
                'Wbar':       <Quantity 0. km / s>,
                'vlsr_error': <Quantity 5. km / s>,
                'max_iter':   1000,
              }
        """
        self.params = params

        self.Ro = params['Ro']
        self.To = params['To']
        self.dTdr = params['dT/dr']
        self.err_v_lsr = params['vlsr_error']

        self.Vo = params['Vo']
        self.Uo = params['Uo']
        self.Wo = params['Wo']

        self.Us = params['Ubar']
        self.Vs = params['Vbar']
        self.Ws = params['Wbar']

        self.max_iter = params['max_iter']

    def calc(self, coords: SkyCoord, v_lsr, farnear):
        """

        :param coords:  astropy SkyCoord object, contains the coordinates of the sources
        :param v_lsr:   astropy Quantity array of v_lsr
        :param farnear: ndarray of farnear
        :return: v_lsr_rev, Dk_best, d_err_high, d_err_low
        """

        v_lsr_rev, Dk_best = self._calc_Dk(coords, v_lsr, farnear, params)

        # params in use
        Ro = self.Ro
        To = self.To
        dTdr = self.dTdr
        err_v_lsr = self.err_v_lsr

        sin_l = np.sin(coords.l.radian)
        R_tp = Ro * sin_l
        T_tp = To + dTdr * (R_tp - Ro)

        del_v = np.zeros_like(farnear) * u.Unit("km/s")
        v_tp = -9999.9 * np.ones_like(farnear) * u.Unit("km/s")

        # quadrant 1
        idx = coords.l.deg < 90
        v_tp[idx] = +T_tp[idx] - To * sin_l[idx]
        idx = np.logical_and(idx, v_lsr_rev > v_tp)
        del_v[idx] = v_lsr_rev[idx] - v_tp[idx]

        # quadrant 4
        idx = coords.l.deg > 270
        v_tp[idx] = -T_tp[idx] - To * sin_l[idx]
        idx = np.logical_and(idx, v_lsr_rev < v_tp)
        del_v[idx] = v_lsr_rev[idx] - v_tp[idx]

        # Adjust v_lsr so as not to exceed tangent point
        v_lsr_adj = v_lsr - del_v

        # get lower velocity distance uncertainty
        v_lsr_low = v_lsr_adj - err_v_lsr

        v_lsr_rev_low, D1 = self._calc_Dk(coords, v_lsr_low, farnear, params)

        idx = abs(v_lsr_rev_low - v_tp) < err_v_lsr
        D1[idx] = Dk_best[idx]

        # get higher velocity distance uncertainty
        v_lsr_high = v_lsr_adj + err_v_lsr
        v_lsr_rev_high, D2 = self._calc_Dk(coords, v_lsr_high, farnear, params)

        idx = abs(v_lsr_rev_high - v_tp) < err_v_lsr
        D2[idx] = Dk_best[idx]

        # calculate + (high) and - (low) errors
        d_err_low = np.empty_like(Dk_best)
        d_err_high = np.empty_like(Dk_best)
        idx = D2 > D1
        d_err_low[idx] = D1[idx] - Dk_best[idx]
        d_err_high[idx] = D2[idx] - Dk_best[idx]

        idx = np.logical_not(idx)
        d_err_low = D2[idx] - Dk_best[idx]
        d_err_high = D1[idx] - Dk_best[idx]

        # If flagged distance (error=0), use other error estimate
        idx = abs(d_err_low) == 0
        d_err_low[idx] = -d_err_high[idx]

        idx = abs(d_err_high) == 0
        d_err_high[idx] = -d_err_low[idx]

        # Check for pathalogical cases...
        idx = Dk_best <= 0
        Dk_best[idx] = 0
        d_err_low[idx] = 0
        d_err_high[idx] = 0

        idx = (Dk_best + d_err_low) <= 0
        d_err_low[idx] = -Dk_best[idx]

        return v_lsr_rev, Dk_best, d_err_high, d_err_low

    def _calc_Dk(self, coords: SkyCoord, v_lsr, farnear, params):

        # params in use
        Vo = self.Vo
        Uo = self.Uo
        Wo = self.Wo
        Ro = self.Ro
        To = self.To
        dTdr = self.dTdr

        Us = self.Us
        Vs = self.Vs
        Ws = self.Ws

        max_iter = self.max_iter

        # get source galactic coordinates
        gal_long_rad = coords.l.radian
        cos_l = np.cos(gal_long_rad)
        sin_l = np.sin(gal_long_rad)

        gal_lat_rad = coords.b.radian
        cos_b = np.cos(gal_lat_rad)
        sin_b = np.sin(gal_lat_rad)

        # convert to true Heliocentric frame (v_rad)
        Uo_IAU = 10.27 * u.Unit('km/s')
        Vo_IAU = 15.32 * u.Unit('km/s')
        Wo_IAU =  7.74 * u.Unit('km/s')

        v_helio = v_lsr - (Vo_IAU * sin_l + Uo_IAU * cos_l) * cos_b - Wo_IAU * sin_b

        # Make "new" V(LSR) using best Solar Motion

        v_newlsr = v_helio + (Vo * sin_l + Uo * cos_l) * cos_b + Wo * sin_b

        # iteration
        del_d = 99 * u.Unit('kpc')
        tol = 0.01 * u.Unit('kpc')
        n_iter = 0
        Dk = 3 * np.ones_like(farnear) * u.Unit('kpc')

        while del_d > tol and n_iter < max_iter:

            Dk_old = Dk

            d_proj = Dk * cos_b

            r_sq = Ro ** 2 + d_proj ** 2 - 2 * Ro * d_proj * cos_l

            r_proj = np.sqrt(r_sq)


            sin_beta = d_proj * sin_l / r_proj
            cos_beta = (Ro - d_proj * cos_l) / r_proj
            beta = np.arctan2(sin_beta, cos_beta)
            beta_deg = beta.to('deg')

            gamma = np.pi - gal_long_rad - beta.value
            cos_gamma = np.cos(gamma)
            sin_gamma = np.sin(gamma)

            v_fixed = v_newlsr - (Vs * sin_gamma - Us * cos_gamma) * cos_b - Ws * sin_b

            Rs = r_proj
            V_proj = v_fixed * cos_b
            D_near, D_far = self._kinematic_distance(V_proj, coords.l, Ro, To, r_proj, dTdr)
            Dk = D_near
            Dk[farnear != 0] = D_far[farnear != 0]

            idx = np.logical_and(D_near.value <= 0, D_far.value > 0)
            Dk[idx] = D_far[idx]

            idx = np.logical_and(D_far <= 0, D_near > 0)
            Dk[idx] = D_near[idx]

            del_d = max(abs(Dk - Dk_old))
            n_iter += 1

        v_lsr_rev = v_fixed

        return v_lsr_rev, Dk

    @staticmethod
    def _kinematic_distance(Vlsr, gal_long, Ro, To, Rs, dTdr):

        glongrad = gal_long.rad

        cos_l = np.cos(glongrad)
        sin_l = np.sin(glongrad)

        Rosinl = Ro * sin_l
        Rocosl = Ro * cos_l

        Theta = To + dTdr * (Rs - Ro)
        Tosinl = To * sin_l
        Tsinl = Theta * sin_l

        rootterm = Rocosl ** 2 + (Tsinl / (Tosinl / Ro + Vlsr / Ro)) ** 2 - Ro ** 2

        rootterm.value[rootterm.value < 0] = 0

        D_near = np.empty_like(gal_long.value) * u.Unit('kpc')
        D_far = np.empty_like(gal_long.value) * u.Unit('kpc')

        for i, gal_l in enumerate(gal_long.deg):

            if  0 <= gal_l < 90:
                D_near[i] = Rocosl[i]- np.sqrt(rootterm[i])
                D_far[i] = Rocosl[i] + np.sqrt(rootterm[i])
            elif 90 <= gal_l < 270:
                D_near[i] = Rocosl[i] + np.sqrt(rootterm[i])
                D_far[i] = D_near[i]
            elif 270 <= gal_l < 360:
                D_near[i] = Rocosl.value[i] - np.sqrt(rootterm[i])
                D_far[i] = Rocosl.value[i] + np.sqrt(rootterm[i])

        return D_near, D_far

    def print_params(self, comment_symbol='!'):
        print("{0} Assumed Galactic parameters: Ro, To, dTdr   {1}, {2}, {3}".format(
               comment_symbol,                       self.Ro, self.To, self.dTdr
        ))
        print("{0} Assumed Solar Motion params: Uo, Vo, Wo     {1}, {2}, {3}".format(
            comment_symbol,                       self.Uo, self.Vo, self.Wo
        ))
        print("{0} Assumed Source Motion params:Us, Vs, Ws     {1}, {2}, {3}".format(
            comment_symbol,                       self.Us, self.Vs, self.Ws
        ))
        print("{0} Assumed uncertainties in V_lsr:             {1}".format(
            comment_symbol,                               self.err_v_lsr
        ))


def dk_output(src_list, coords, v_lsr, v_lsr_rev, Dk_best, d_err_high, d_err_low, output_filename=None, comment_symbol: str='!'):
    output_df = pd.DataFrame()
    output_df['Source'] = src_list
    output_df['Gal Long'] = coords.l.deg
    output_df['Gal Lat'] = coords.b.deg
    output_df['V_lsr'] = v_lsr
    output_df['V_rev'] = v_lsr_rev
    output_df['Rev. D_k'] = Dk_best
    output_df['+'] = d_err_high
    output_df['-'] = d_err_low

    formatters = [
        lambda x: "{:12}".format(x),
        lambda x: "  {:.3f}".format(x),
        lambda x: "  {:.3f}".format(x),
        lambda x: "    {:.1f}".format(x),
        lambda x: "    {:.1f}".format(x),
        lambda x: "    {:.2f}".format(x),
        lambda x: "  {:.2f}".format(x),
        lambda x: " {:.2f}".format(x),
    ]

    headers = [
        'Source       ',
        'Gal Long  ',
        'Gal Lat    ',
        'V_lsr     ',
        'V_rev    ',
        'Rev. D_k     ',
        '+/-'
    ]

    units = [
        comment_symbol,
        ' ' * (15 - len(comment_symbol)),
        '  (deg)',
        '    (deg)',
        '    (km/s)',
        '    (km/s)',
        '     (kpc)',
        '      (kpc)',
        ]
    if output_filename is None:
        print(comment_symbol, ''.join(headers))
        print(''.join(units))
        print(output_df.to_string(header=False, index=False, formatters=formatters, max_rows=None))
    else:
        output_df.to_csv(args.output, sep=',', float_format="%0.5f", index=False)


if __name__ == '__main__':
    args = parse_args()

    params = get_parameters_from_json(args.parameter_file)
    params['max_iter'] = args.max_iter
    comment_symbol = args.comment_symbol
    src_list, coords, v_lsr, farnear = get_source_data_from_file(args.source_file,
                                                                 args.use_radec,
                                                                 data_file_comment_symbol=comment_symbol)

    dk_calculator = DkCalculator(params)

    v_lsr_rev, Dk_best, d_err_high, d_err_low = dk_calculator.calc(coords, v_lsr, farnear)

    dk_calculator.print_params(comment_symbol)

    dk_output(src_list,
              coords,
              v_lsr,
              v_lsr_rev,
              Dk_best,
              d_err_high,
              d_err_low,
              args.output,
              comment_symbol=comment_symbol)

