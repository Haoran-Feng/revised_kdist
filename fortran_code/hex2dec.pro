; convert hexical format to decimal format
function hex2dec,input,hour=hour,reverse=reverse,mksign=mksign,decimal=decimal,separator=separator

;usage hex2dec x=hex2dec(input[,/hour])
if(n_params() lt 1) then begin
print,"usage hex2dec x=hex2dec(input[,/hour,/reverse])"
return,0
endif

if n_elements(separator) then sp=separator else sp=':'

if (not keyword_set(reverse)) then begin
;strl=strlen(input)
;smark=strmid(input,0,1)
;if(smark[0] eq "+") then begin
;	sign=1
;	inpstr=strmid(input,1,strl-1)
;endif else begin
;	if (smark[0] eq "-") then begin
;	sign=-1
;	inpstr=strmid(input,1,strl-1)
;	endif else begin
;	sign=1
;	inpstr=input
;	endelse
;endelse

strl=strlen(input)

;col1=stregex(inpstr,':')
;col2=stregex(strmid(inpstr,col1+1,strl-col1-1),':')+col1+1
ssa=strsplit(input,sp,/extract)
rst=[]
;print,n_elements(input)
for i1=0,n_elements(input)-1 do begin
    if isa(ssa,'list') then ssb=ssa[i1] else ssb=ssa
    
for i=0,n_elements(ssb)-1 do begin
    ss=ssb

if ss[0] eq '+' or ss[0] eq '-' then begin
    tmp=[strcompress(ss[0]+ss[1],/remove),ss[2:*]]
    ss=tmp
    endif
narea=n_elements(ss)
ss0=ss[0]
ss=float(ss)
if strmid(ss0,0,1) eq '-' then sg=-1 else sg=1
;if ss[0] ge 0 then sg=1 else sg=-1
ss=abs(ss)
;print,ss
;hdeg=double(strmid(inpstr,0,col1))
;if (col2 gt col1) then begin
;hm=double(strmid(inpstr,col1+1,col2-col1-1))
;hs=double(strmid(inpstr,col2+1,strl-col2-1))
;endif else begin
;hm=double(strmid(inpstr,col1+1,strl-col1-1))
;hs=0.0D
;endelse

;print,hdeg,hm,hs
x=0
for i=0,narea-1 do $
x=(x+ss[narea-1-i])/60.0
x=x*sg*60
if(keyword_set(hour)) then x=x*15.0
rst=[rst,x]
endfor
endfor
;print,x
return, rst

endif else begin
rst=[]
;print,n_elements(input)
for i1=0,n_elements(input)-1 do begin  
ss=input[i1]      
if((size(ss,/tname) ne 'DOUBLE') and (size(ss,/tname) ne 'FLOAT') ) then begin
print,'input data type error'
return,0
endif
if (keyword_set(hour)) then opt=ss/15.0 else opt=ss
if(ss ge 0.0) then sign='+' else sign='-'
opt=abs(opt)
dd=fix(opt)
mm=fix((opt-dd)*60)
ss=((opt-dd)*60-mm)*60

dds=string(dd,format='(I2)')
if (dd lt 10.0) then dds='0'+strmid(dds,1,1)
mms=string(mm,format='(I02)')
;if (mm lt 10.0) then mms='0'+strmid(mms,1,1)
if n_elements(decimal) then fmt=strcompress('(F0'+string(decimal+3)+'.'+string(decimal)+')',/remove) else fmt='(I02)'
sss=string(ss,format=fmt)
;if (ss lt 10.0) then sss='0'+strmid(sss,1,4)
if(keyword_set(mksign) or (sign eq '-')) then ostr=sign+dds+sp+mms+sp+sss else $
ostr=dds+sp+mms+sp+sss

rst=[rst,ostr]
endfor
return,rst

endelse

end