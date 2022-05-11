PRO read_bon_fullfile,bonn_filename,azm0,range0,elev0,x0,y0,h0,zh0,zdr0,rho0,phi0,kdp0,vr0,ngates1,tstart

ngates=933
nbeams=360
nelev=10
zh0=fltarr(nelev,ngates,nbeams)
zdr0=fltarr(nelev,ngates,nbeams)
rho0=fltarr(nelev,ngates,nbeams)
phi0=fltarr(nelev,ngates,nbeams)
kdp0=fltarr(nelev,ngates,nbeams)
vr0=fltarr(nelev,ngates,nbeams)
h0=fltarr(nelev,ngates,nbeams)
x0=fltarr(nelev,ngates,nbeams)
y0=fltarr(nelev,ngates,nbeams)
bonn_elevation=fltarr(nelev)
bonn_filename0=strarr(nelev)
azm0=fltarr(nelev,nbeams)
range0=fltarr(nelev,ngates)
tstart=strarr(nelev)
tend=strarr(nelev)
ngates1=fltarr(nelev)

date=file_basename(bonn_filename,'_00.h5')
print,date
datentime=strmid(date,19,12)
aux=strmid(datentime,8,4)
aux1=strmid(datentime,8,3)
time=float(aux)
date=strmid(datentime,0,8)

; Which moments should be read?
moments=['Zh','ZDR','Vh','RHOHV','PHIDP','KDP','UZh']

bonn_elevation(0)=1.
bonn_elevation(1)=2.
bonn_elevation(2)=3.1
bonn_elevation(3)=4.5
bonn_elevation(4)=6.0
bonn_elevation(5)=8.2
bonn_elevation(6)=11.
bonn_elevation(7)=14.
bonn_elevation(8)=18.
bonn_elevation(9)=28.

filename=file_search('/automount/radar/scans/2017/2017-07/2017-07-06/n_ppi_*/n_ppi_*deg_12345_'+date+'*_00.h5')
;filename=file_search('/automount/radar/scans/2015/2015-01/2015-01-14/n_ppi_*/*0*.mvol')

files=file_basename(filename,'_00.h5')
datentimefiles=strmid(files,19,12)
times=strmid(datentimefiles,8,4)
delta_t=abs(times-time(0))

elev=strmid(filename,47,9)

aux2=min(delta_t(where(elev eq 'n_ppi_010')),index)
if times(index)-time(0) lt 0 then index=index+1
bonn_filename0(0)=filename(index)
index0=where(elev eq 'n_ppi_020')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(1)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_031')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(2)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_045')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(3)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_060')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(4)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_082')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(5)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_110')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(6)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_140')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(7)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_180')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(8)=filename(index+index0(0))
index0=where(elev eq 'n_ppi_280')
aux2=min(delta_t(index0),index)
if times(index+index0(0))-time(0) lt 0 then index=index+1
bonn_filename0(9)=filename(index+index0(0))

filetime=fltarr(nelev)
for ielev=0,nelev-1 do begin

read_boxpol_aziscan_new, filename=bonn_filename0(ielev), moments=moments,data=bonn_data,azi=bonn_azi,elevation=bonn_elevation(ielev),range=bonn_range,time=bonn_time,verbose=1,turn_scan=turn_scan,offset_range=150.,prf1=bonn_prf1,prf2=bonn_prf2,tsamples=bonn_tsamples,scanspeed=bonn_scanspeed,date_start=date_start,date_end=date_end;

zzz=bonn_data(*,*,0)
zdr=bonn_data(*,*,1)
vr=bonn_data(*,*,2)
rho=bonn_data(*,*,3)
phi=bonn_data(*,*,4)

aux=file_basename(bonn_filename0(ielev),'_00.h5')
filetime(ielev)=strmid(aux,27,6)

radar_height=99.5 ;in m
radar_height=radar_height/1000.
lonlat_rad0=[7.071663,50.73052]
lat_rad=lonlat_rad0(1)
lon_rad=lonlat_rad0(0)
azm=bonn_azi
distance=bonn_range/1000.

ngates0=n_elements(bonn_range)
nbeams0=n_elements(bonn_azi)

ngates1(ielev)=ngates0
enhance_phidp, data=bonn_data,moments=moments,ephidp=phi,kdp=kdp,calculate_kdp=1,fill_gaps=0,range=bonn_range,phi_gaps=phi_gaps,smoothekm=3000,calculate_phi0=1

;*******************************************************************
ephidp=phi

; Correction earth radius
Rcorr=6370.*4/3
s=fltarr(ngates0,nbeams0)
s0=fltarr(ngates0,nbeams0)
h=fltarr(ngates0,nbeams0)
snr=fltarr(ngates0,nbeams0)
lonlat=fltarr(ngates0,nbeams0,2)
x=fltarr(ngates0,nbeams0)
y=fltarr(ngates0,nbeams0)


for i=0,ngates0-1 do begin
    for j=0,nbeams0-1 do begin
        s(i,j)=Rcorr*atan(distance(i)*cos(bonn_elevation(ielev)*!pi/180.)/(Rcorr+distance(i)*sin(bonn_elevation(ielev)*!pi/180.)))
        h(i,j)=sqrt(distance(i)^2+2*distance(i)*Rcorr*sin(bonn_elevation(ielev)*!pi/180.)+Rcorr^2)-Rcorr+radar_height
        x(i,j)=s(i,j)*sin(azm(j)*!pi/180.)
        y(i,j)=s(i,j)*cos(azm(j)*!pi/180.)
    endfor
endfor

; attenuation correction for surface Zh and Zdr
; to do attenuation correction, Ryzhkov et al., 2013: Polarimetric Radar
; characteristics of melting hail. Part II
; code by Xinxin Xie, originally in matlab

alpha = 0.3; % dB/deg
b = 0.5;
a = 2.6e-3; 
beta = 0.03

zh = zzz;
zdr = zdr;
zh_bias = 0.5;
zdr_bias = 0.8;
zh1 = 10.^(0.1*(zh));
Ahh = zh*0;
delta_zdr = zdr*0;
zh_corr=zh;fltarr(ngates,nbeams)
zdr_corr=zdr;fltarr(ngates,nbeams)
rm = 750-1;
if bonn_elevation(ielev) eq 6. then rm=932.
if bonn_elevation(ielev) eq 8.2 then rm=799.
if bonn_elevation(ielev) eq 11 then rm=639.
if bonn_elevation(ielev) eq 14 then rm=495.
if bonn_elevation(ielev) eq 18 then rm=399.
if bonn_elevation(ielev) eq 28 then rm=287.
resolution=bonn_range(2)-bonn_range(1)
print, resolution

for i = 0,nbeams0-1 do begin
    zh1(0:rm,i) = 10.^(0.1*(zh(0:rm,i)+zh_bias));

    row = where(zh(0:rm,i) gt 50);
    if n_elements(row) gt 2 then begin
        r1 = row(0);
        r2 = row(n_elements(row)-1);
        rowj = 0;
        for r = 1,n_elements(row)-1 do begin
            rowx = row(r)-row(r-1);
            if (rowx gt 5  and rowj eq 0) then begin
                r2 = row(r-1);
                rowj = 1;
            endif
        endfor

        for r = r1,r2 do begin
            I_r1r2 = 0.46*b*total(zh1(r1:r2,i)^b*resolution*1e-3);
            I_rr2 =  0.46*b*total(zh1(r:r2,i)^b*resolution*1e-3);
            PIA = 2*total(a*zh1(r1:r2,i)^b*resolution*1e-3);
            C = exp(0.23*PIA*b)-1;
            Ah = zh1(r,i)^b* C / (I_r1r2 + C*I_rr2); 
            Ahh(r,i) = Ah; 
            zh_corr(r,i) = zh(r,i) +zh_bias+ 2*total(Ahh(r1:r,i)*resolution*1e-3);
        endfor
        for r = r2,rm do begin
            delta_z = 0.28*(phi(r,i)-phi(r2,i));
            zh_corr(r,i) =  zh(r,i) + zh_bias + delta_z;
        endfor
        row = where(zh_corr(0:rm,i) gt 50);
        r1 = row(0);
        r2 = row(n_elements(row)-1);
        rowj = 0;
        for r = 1,n_elements(row)-1 do begin
            rowx = row(r)-row(r-1);
            if (rowx gt 5 and  rowj eq 0) then begin
                r2 = row(r-1);
                rowj = 1;
            endif
        endfor

        row1 = where( zh_corr(r2:rm,i) lt 30 and  zh_corr(r2:rm,i) gt 20 )+r2-1;
        if n_elements(row1) ge 1 then begin
            r3 = row1(0);
            r4 = row1(n_elements(row1)-1);
            rowj = 0;
            for r = 1,n_elements(row1)-1 do begin
                rowx = row1(r)-row1(r-1);
                if (rowx gt 10 and rowj eq 0) then begin
                   r4 = row1(r-1);
                   rowj = 1;
                endif
            endfor

            for r = r3,r4  do begin
                zdr_corr(r,i) = 0.585 - 0.0507*zh_corr(r,i)+0.00165*zh_corr(r,i)*zh_corr(r,i);  
            endfor
            delta_zdr1 = mean(abs(zdr(r3:r4,i) + zdr_bias - zdr_corr(r3:r4,i)));
;            beta_h = max(delta_zdr1 / (phi(r2,i) - phi(r1,i)));
            beta_h = 0.11
            for r = 0,r1 do begin
                if zh_corr(r,i) gt 0 then begin
                   delta_zdr(r,i) = beta*(phi(r)-phi(0));
                   zdr_corr(r,i) = zdr(r,i) + zdr_bias + delta_zdr(r,i);
                endif
            end
            for r = r1,r2 do begin
                if zh_corr(r,i) gt 0 then begin
                   delta_zdr(r,i) = beta_h*(phi(r,i)-phi(0,i)) - (beta_h-beta)*(phi(r1,i)-phi(0,i));
                   zdr_corr(r,i) = zdr(r,i) + zdr_bias + delta_zdr(r,i); 
                endif
            endfor
            for r = r4,rm  do begin
                zdr_corr(r,i) = zdr(r,i) + zdr_bias;
                if rho(r,i) gt 0.8 and zh_corr(r,i) gt 0 then begin
                   delta_zdr(r,i) = beta*(phi(r,i)-phi(0,i)) + (beta_h-beta)*(phi(r2,i)-phi(r1,i));
                   zdr_corr(r,i) = zdr(r,i)+zdr_bias+delta_zdr(r,i); 
                endif
            endfor
            for r = r2,r3 do begin
                if zh_corr(r,i) gt 0 then begin
                    delta_zdr(r,i) = beta*(phi(r,i)-phi(0,i)) + (beta_h-beta)*(phi(r2,i)-phi(r1,i));
                    zdr_corr(r,i) = zdr(r,i)+zdr_bias+delta_zdr(r,i); 
                endif
            endfor

        endif
    endif else begin
        zh_corr(*,i)=zh(*,i)+0.28*(phi(*,i))+zh_bias
        zdr_corr(*,i)=zdr(*,i)+0.03*(phi(*,i))+zdr_bias
    endelse
endfor

h0(ielev,0:ngates0-1,*)=h(*,*)
x0(ielev,0:ngates0-1,*)=x(*,*)
y0(ielev,0:ngates0-1,*)=y(*,*)
zh0(ielev,0:ngates0-1,*)=zh_corr(*,*)
zdr0(ielev,0:ngates0-1,*)=zdr_corr(*,*)
rho0(ielev,0:ngates0-1,*)=rho(*,*)
phi0(ielev,0:ngates0-1,*)=phi(*,*)
kdp0(ielev,0:ngates0-1,*)=kdp(*,*)
vr0(ielev,0:ngates0-1,*)=vr(*,*)
azm0(ielev,*)=azm(*)
range0(ielev,0:ngates0-1)=distance(*)
elev0=bonn_elevation
tstart(ielev)=date_start
tend(ielev)=date_end

endfor		;nelev

end



