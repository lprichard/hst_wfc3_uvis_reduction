; program to figure out which darks are needed for a set of
; observations. Useful to make sure you process the darks before you
; need them

pro find_darks, postflash=postflash, smooth=smooth, gzip=gzip, avg=avg
if not keyword_set(gzip) then  spawn, 'ls -1 *raw.fits',raw else spawn, 'ls -1 *raw.fits.gz',raw

  len=n_elements(raw)
  darkfilearr = strarr(n_elements(raw))
  filterlist = strarr(n_elements(raw))
  rootarr =  strarr(n_elements(raw))
  detectorlist = strarr(n_elements(raw))

  darklookupfilearr = strarr(n_elements(raw))

  for i=0,len-1 do begin
     name=raw(i)
     h=headfits(name) 
     filter=strcompress(sxpar(h,'FILTER'),/remove_all)    
     detector=strcompress(sxpar(h,'DETECTOR'),/remove_all)  
     dateobs = strcompress(sxpar(h,'DATE-OBS'),/remove_all)  
     timeobs = strcompress(sxpar(h,'TIME-OBS'),/remove_all)  
     useafter=strcompress(sxpar(h,'DETECTOR'),/remove_all)  
     ; modified julian date exposure end time
     EXPEND = double(strtrim(sxpar(h, 'EXPEND'),2))

     newdark = read_dark_lookup(expend,smooth=smooth,postflash=postflash, avg=avg)

     detectorlist[i] = detector
     filterlist[i] = filter

     if detector eq 'UVIS' then begin
        dark=strmid(strcompress(sxpar(h,'DARKFILE'),/remove_all),5,13)
        darkfilearr[i] = dark
        darklookupfilearr[i] = newdark

        root= strmid(name, 17, 9, /reverse_offset)
        rootarr[i] = root
        
        print, root, '   ', filter, '    ', dark, '    ', dateobs, '    ', timeobs, '    ', newdark
        
     endif
     
  endfor
     ; Lets figure out which darks there are
; need to sort first, then use the uniq function on sorted array. That
; gives you the subscripts. 
        darkfiles = darkfilearr[sort(darkfilearr)]
        diffdark = uniq(darkfiles)
        
        darksneeded = darkfiles[diffdark]
        print, '---------------------------------------------------------'
        print, 'Make sure you process the following darks:'
        print, darksneeded
        print, '---------------------------------------------------------'
        ldarkfiles = darklookupfilearr[sort(darklookupfilearr)]
        difflookupdark = uniq(ldarkfiles)
        
        darkslookupneeded = ldarkfiles[difflookupdark]


        print, '---------------------------------------------------------'
        print, 'These will be called the following (based on /smooth, and /postflash, and dark lookup table):'
        print, darkslookupneeded
        print, '---------------------------------------------------------'


  
end
