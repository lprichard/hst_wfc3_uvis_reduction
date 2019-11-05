; program to read dark_lookup tables for UVIS WFC3
; Takes a Julian date, and returns filename for dark for given julian date.
; binned option calls the binned files
; /avg uses the averaged darks after the smooth iteration

function read_dark_lookup, julend, binned=binned, smooth=smooth, postflash=postflash, avg=avg, dkc=dkc

  dir = '~/candels/calibrations/'
;  if keyword_set(binned) then begin
;     if keyword_set(smooth) then file = 'dark_lookup_binned_smoothed.txt'  else file = 'dark_lookup_binned.txt' 
;     endif else begin
;        if keyword_set(smooth) then begin
;           file='dark_lookup_unbinned_smoothed.txt' 
;           ; by default, postflash is not also binned
;           if keyword_set(postflash) then file='dark_lookup_postflash_smoothed.txt' 
;        endif else begin
;           file='dark_lookup_unbinned.txt'
;        endelse
;     endelse

  if keyword_set(dkc) then file = 'dark_lookup_dkc.txt' else file='dark_lookup.txt'
 
  thisformat = '(A19, A3, I2, I4, A8)'

  readcol, dir+file, name, mon, day, yr, time, format=thisformat, /silent
  
  for idx=0, n_elements(name)-1 do begin
     temp = name[idx]
     if keyword_set(smooth) then strput, name, 's', 1
     if keyword_set(avg) then strput, name, 'a', 1
     if keyword_set(binned) then strput, name, 'b', 0
     if keyword_set(postflash) then strput, name, 'p', 0
     name[idx] = temp
  endfor

  juldate = dblarr(n_elements(name))
  for idx=0, n_elements(name)-1 do begin
   endswith = strmid(name[idx],7, /reverse_offset)

  sday =strtrim(string(day[idx]),2)
  syr = strtrim(string(yr[idx]),2)

     date = sday+'-'+mon[idx]+'-'+syr+' '+ time[idx]
     juldate[idx] = date_conv(date, 'M')
     ;print, date
     ;print, juldate[idx]
  endfor
  
; This finds the closest dark AFTER the use after date.
diff = julend-juldate
neg = where(diff le 0)
diff[neg] = 999.99
minval = min(diff, use)

return, name[use]

end
