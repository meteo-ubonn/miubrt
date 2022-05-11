pro enhance_phidp, data=data,moments=moments,ephidp=phidp,kdp=kdp,calculate_kdp=calculate_kdp,fill_gaps=fill_gaps,range=range,phi_gaps=phi_gaps,phi0=phi0,smoothekm=smoothekm,calculate_phi0=calculate_phi0,max_phi0range=max_phi0range,totalz=totalz,calculated_phi0=calculated_phi0,phi0azi=phi0azi

; This subroutine creates a filtered and smoothed set of differential phase PHIdp estimates. If calculate_kdp is set to one, kdp is also returned

; Data: Must contain the polarimetric radar measurements of the parameters ZH,UH,PHIdip,RHOhv  in the form   data(range,azimuth,parameter)
; moments: must contain an array of strings indicating which parameter index of data corresponds to which parameter
; ephidp : returns the filtered and smoothed estimates of PHIdp
; kdp: returns estimates of kdp if calculate_kdp is set to 1
; fill_gaps: if set to 1, the places where no PHIdp was measured are filled the last measurement PHIdp (in directon of tha radar)
; range: if calculate_kdp is set to one, this mut contain the ranges from the radar (to convert to degrees per kilometer). Not needed if kdp is not calcaulated
; phi_gaps= also returns the filtered estimates of phidp, but if "fill_gaps" is set to one, this phidp is still returned with gaps.

;!except=0

print,'this is the routine being used here'
if n_elements(calculate_kdp) eq 0 then calculate_kdp=0
if n_elements(calculate_phi0) eq 0 then calculate_phi0=0
if n_elements(fill_gaps) eq 0 then fill_gaps=0

nrange=n_elements(data(*,0,0))
na=n_elements(data(0,*,0))

izh=where(moments eq 'zh')
iuh=where(moments eq 'uh')
iphidp=where(moments eq 'phidp')
irhohv=where(moments eq 'rhohv')
izdr=where(moments eq 'zdr')

;izh=where(moments eq 'ZH')
;iuh=where(moments eq 'UZh')
;iphidp=where(moments eq 'PHIDP')
;irhohv=where(moments eq 'RHOHV')
;izdr=where(moments eq 'ZDR')

kdps=data(*,*,iphidp)
kdps(1:nrange-1,*)=kdps(1:nrange-1,*)-kdps(0:nrange-2,*)

kdps2=data(*,*,iphidp)
kdps2(0:nrange-2,*)=kdps2(1:nrange-1,*)-kdps2(0:nrange-2,*)


kickphidp=where(data(*,*,iuh)-data(*,*,izh) gt 2 or data(*,*,irhohv) lt 0.87 or kdps lt -15 or kdps gt 15 or kdps2 lt -15 or kdps2 gt 15 or data(*,*,izdr) gt 5.5 or data(*,*,izdr) lt -6 or data(*,*,izdr) gt 1.42+0.0667*data(*,*,izh)+0.000485*data(*,*,izh)^2.  or data(*,*,izh) lt 5)
;kickphidp=where(data(*,*,irhohv) lt 0.8)

phidp=data(*,*,iphidp)
if kickphidp(0) ne -1 then phidp(kickphidp)=0./0. 

;phidp=smooth(phidp,3,/nan)

pval=fltarr(nrange,na)
nval=fltarr(nrange,na)

len=(smoothekm/(range(1)-range(0)))/2.
len2=len/1.5

print, smoothekm,len, len2,range(1)-range(0)

phi0azi=fltarr(na)
n=fltarr(na)
phidp2=phidp
nfin=fltarr(nrange,na)
for ai=0,na-1 do begin
for ri=0,nrange-1 do begin
  ;if finite(phidp(ri,ai)) then begin
  up=(ri+len < nrange-1)
  down=(ri-len > 0)
  fin=where(finite(phidp2(down:up,ai)),nf)
  if nf gt 0 then begin nfin(ri,ai)=total(data(down+fin,ai,irhohv))
   ;if finite(phidp2(ri,ai)) and nfin(ri,ai) gt len2  then phidp(ri,ai)=median(phidp2(down:up,ai))
   if nfin(ri,ai) gt len2 then phidp(ri,ai)=median(phidp2(fin+down,ai))
  endif
endfor
phidp(*,ai)=smooth(phidp(*,ai),3,/nan)

for ri=0,nrange-1 do begin
  if n(ai) lt 6 and finite(phidp(ri,ai)) and phidp(ri,ai)  lt -50 then begin
    n(ai)=n(ai)+1
    phi0azi(ai)= phi0azi(ai)+phidp(ri,ai)
  endif

endfor

endfor
phi0azi=phi0azi/n


  ngood=0
  ri=3
  while ngood lt 500  and ri lt nrange-1 do begin
    good=where( finite(data(2:ri,*,iphidp)) and data(2:ri,*,izdr) lt 3 and data(2:ri,*,izdr) gt -3 and data(2:ri,*,irhohv) gt 0.97,ngood   )
    ri=ri+1
  endwhile
  testvals=data(2:ri-1,*,iphidp)
  max_phi0range=range(ri-1)
  dummyz=data(2:ri-1,*,izh)
  if ngood ne 0 then totalz=total(dummyz(good)) else totalz=-32.
  ;phi0=total(testvals(good))/ngood
  ;inliers=where( (testvals-phi0)^2 lt 20.^2.,ngood)
  ;phi0=total(testvals(inliers))/ngood
  nogood=0
	if ngood ne 0 then calculated_phi0=median(testvals(good)) else begin
	     calculated_phi0=-92
	     nogood=1
	     print, 'NO PHI0 FOUND, SETTING ALL PHIDP TO 0!!!!!!!!!!!!!!!!!!!'
	     
	 endelse    
print, 'calculated phi0 ', calculated_phi0,ri 

if calculate_phi0 then phi0=calculated_phi0

phidp=phidp-phi0

for ai=0,na-1 do begin
for ri=0,10 do if phidp(ri,ai) gt 10 then phidp(ri,ai)=0./0.
for ri=10,20 do if phidp(ri,ai) gt 20 then phidp(ri,ai)=0./0.

lphi=0
lconf=len2
lra=0
for ri=00,nrange-1 do begin
  if finite(phidp(ri,ai)) then begin
    if phidp(ri,ai) lt lphi-1 or phidp(ri,ai) gt lphi+5 then begin
     if nfin(ri,ai) gt lconf then begin
      lphi=phidp(ri,ai)
      lconf=nfin(ri,ai)
      phidp(lra,ai)=0./0.
     endif else phidp(ri,ai)=0./0.
    endif else begin
      lphi=phidp(ri,ai)
      lconf=nfin(ri,ai)   
      lra=ri
    endelse
  endif
endfor

endfor

;for ai=0,na-1 do if finite(phi0azi(ai)) and abs(phi0azi(ai)-phi0) lt 10  then phidp(*,ai)=phidp(*,ai)-phi0azi(ai)+phi0
;if n_elements(
;restore, file='/Data2/Off_phi0


print, 'phi offset: ',phi0
bad=where(phidp lt 0)
if bad(0) ne -1 then phidp(bad)=0.

;if fill_gaps then begin
  phi_gaps=phidp
  lastfinite=intarr(na)
	for ai=0,na-1 do begin
		for ri=1,nrange-1 do begin
      if finite(phidp(ri,ai)) and phidp(ri,ai) gt 0 then begin
      	if lastfinite(ai) eq 0 then begin
      	  phidp(0:ri-1,ai)=0;phidp(ri,ai)
      	endif else begin
      	  if ri-lastfinite(ai) gt 1 then begin
      	     phidp(lastfinite(ai)+1:ri-1,ai)= (phidp(ri,ai)*nfin(ri,ai) + phidp(lastfinite(ai),ai)*nfin(lastfinite(ai),ai))/( nfin(lastfinite(ai),ai)+nfin(ri,ai)  )
      	        
      	  endif   
      	endelse
      	lastfinite(ai)=ri
      endif else begin
        ;if ri ge 2 then begin
        ;  if phidp(ri-1,ai)-phidp(ri-2,ai) lt 20 then phidp(ri,ai)=phidp(ri-1,ai) else phidp(ri-1:ri,ai)=phidp(ri-2,ai)
        ;endif else begin
        ;  phidp(ri,ai)=phidp(ri-1,ai)
        ;endelse
      endelse  
		endfor
		if lastfinite(ai) ne ri-1 then phidp(lastfinite(ai)+1:ri-1,ai)= phidp(lastfinite(ai),ai)
		phidp(*,ai)=smooth(phidp(*,ai),3,/nan)
		
		if lastfinite(ai) eq 0 then phidp(*,ai)=0.
	endfor
	;phidp(*,ai)=smooth(phidp(*,ai),3,/nan)
	
;endif

if 0 then begin
window,0

for ai=0,na-1 do begin
  oplot, [0,max(range)],[0,0]
  plot, range,data(*,ai,iphidp)-phi0,xrange=[0,max(range)],xstyle=1,yrange=[-180,180],ystyle=1,psym=3

  oplot,range, phidp(*,ai),psym=3,color=180
  oplot, range,data(*,ai,iphidp)-phi0,psym=3
  oplot, range,data(*,ai,izh)-100.,psym=3
   oplot, range,data(*,ai,irhohv)*20.-150. ,color=180,psym=3
  print, ai
  halt=get_kbrd(1)
endfor
endif

if calculate_kdp then begin


bin=(range(1)-range(0))/1000.
 
smootherange=smoothekm/(range(1)-range(0))

s1=intarr(nrange)
s2=intarr(nrange)

for ri=0,nrange-1 do begin
	s1(ri)=max([0,ri-smootherange/2])
	s2(ri)=min([ri+smootherange/2,nrange-1])
endfor		
	
kdp=fltarr(nrange,na)
lastfinite=lastfinite+smootherange/2
lastfinite=lastfinite*(lastfinite lt nrange)+nrange*(lastfinite ge nrange) < nrange-1
for ai=0,na-1 do begin
   for ri=0,lastfinite(ai)-1 do begin
		x=2.*bin*indgen(s2(ri)-s1(ri)+1)
		y=phidp(s1(ri):s2(ri),ai)
		aux=linfit(x,y)
        kdp(ri,ai)=aux(1)
   endfor
endfor
print,'here is the kdp calculation'

endif

if nogood then phidp=phidp*0.

end
