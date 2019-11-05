; Program to list the raw darks that are needed to download for a set
; of darks.

; To use this to list all darks easily:

;find_rawdarks, list='batch2.list', /CDBS, /all
; if you want a single dark, set single to that value.

; newdarks handles the new version 2.0 darks 
; ctedarks handles the new cte corrected version 2.0 darks

; find_rawdarks, /newdarks, /all
; find_rawdarks, list='batch51.list', /all, /cdbs, /new


pro find_rawdarks, list=list, all=all, CDBS=CDBS, single=single, newdarks=newdarks, ctedarks=ctedarks
  
  if not keyword_set(single) then begin
     if keyword_set(list) then begin
        readcol, list, drk, format='A'
     endif else begin
        if not keyword_set(ctedarks) then spawn, 'ls -1 *drk.fits',drk else spawn, 'ls -1 *dkc.fits',drk
     endelse
  endif else drk = single

  len=n_elements(drk)

  for idx=0, len-1 do begin
     
     if not keyword_set(all) then begin
        print, drk[idx]
        print, ''
     endif
     
     if keyword_set(CDBS) then begin
        cdbsdir = '/Users/mrafelski/ThunderBay/rawdarks/CDBS/drk_files/'
        hdr = xheadfits(cdbsdir+drk[idx])
     endif else begin
        hdr = xheadfits(drk[idx])
     endelse
     
     if not keyword_set(newdarks) then fhdr = strpos(hdr, 'HISTORY nn') else fhdr = strpos(hdr, 'HISTORY id')
     
     gd = where(fhdr eq 0)
     
     shdr = hdr[gd]
     
;print, shdr
     
   rawdarks = shdr
   
   
   for idy=0, n_elements(shdr)-1 do begin
      if not keyword_set(newdarks) then rawdarks[idy] = strmid(rawdarks[idy],10,9)
      if keyword_set(newdarks) then rawdarks[idy] = strmid(rawdarks[idy],8,9)

;rawdarks[idy] = strmid(rawdarks[idy],8,24)
      
   endfor

   if keyword_set(newdarks) then begin
      gd = where(rawdarks ne 'identifyi')
      rawdarks=rawdarks[gd]
   endif
   
   forprint, rawdarks
   if not keyword_set(all) then print, ''
   
endfor
  
end
