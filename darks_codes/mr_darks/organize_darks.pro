; If you downloaded a bunch of data files and want to organize them
; based on a list of input files, this is how you do it. 

; This is written without testing. Test it, and first check directory
; structure is right. 

; listofdarks is the filename of the CDBS darks you want to simulate
; inputdatadir is where the data are located.
; ctecor lets you ctecorrect the data with this code as well. If ctecor is not called, then assumes data is already cte corrected

; organize_darks, 'newdarks1.list', /ctecor

; version with no mp possible (e.g. tiger)
; precte means copy files into precte. 
; without it, it assumes cte corrections have been made. 
; removing gzip for now.  Causes problems with darks creation - need
;                                                               to
;                                                               modify code.

; use /precte to 

pro organize_darks, listofdarks, inputdatadir=inputdatadir, precte=precte, newdarks=newdarks

  if not keyword_set(newdarks) then begin
     print, 'Are these not new darks? use /newdarks'
     stop
  endif
;CDBS_dir = '~/edata2/rawdarks/cdbs/drk_files/'
;dataloc = '~/edata2/rawdarks/process/'
;orig_dataloc = '~/edata2/rawdarks/'
;precte_dataloc = '~/edata2/rawdarks/precte/'
CDBS_dir = '~/ThunderBay/rawdarks/cdbs/drk_files/'
dataloc = '~/ThunderBay/rawdarks/process/'
orig_dataloc = '~/ThunderBay/rawdarks/'
precte_dataloc = '~/ThunderBay/rawdarks/precte/'

;if not keyword_set(inputdatadir) then inputdatadir = '~/edata2/rawdarks/download/batch16/'
;if not keyword_set(inputdatadir) then inputdatadir = '/Volumes/ThunderBay/rawdarks/download/batch37_precte/'
if not keyword_set(inputdatadir) then inputdatadir = '/Volumes/ThunderBay/rawdarks/download/batch52/'
;if not keyword_set(inputdatadir) then inputdatadir = '/Users/mrafelski/hststage/batch51_precte/'
;inputdatadir = '/Volumes/LaCie/rawdarks/download/batch17_precte/'
readcol, listofdarks, darkfiles, format='(A)'

cd, inputdatadir, current=origdir

;if not keyword_set(precte) then spawn, 'gunzip *.gz'
;spawn, 'gunzip *.gz'

   darkfiles_root = darkfiles

for idx=0, n_elements(darkfiles)-1 do begin

   boo = strmid(darkfiles[idx], 0, 9)
   darkfiles_root[idx] = boo

   superdark = CDBS_dir+darkfiles[idx]

   hdr = xheadfits(superdark)
   ;fhdr = strpos(hdr, 'HISTORY nn')
   if not keyword_set(newdarks) then fhdr = strpos(hdr, 'HISTORY nn') else fhdr = strpos(hdr, 'HISTORY id')
   gd = where(fhdr eq 0)

   shdr = hdr[gd]
   rawdarks = shdr
   for idy=0, n_elements(shdr)-1 do begin
      if not keyword_set(newdarks) then rawdarks[idy] = strmid(rawdarks[idy],10,9) + '_raw.fits'
      if keyword_set(newdarks) then rawdarks[idy] = strmid(rawdarks[idy],8,9) + '_raw.fits'

      if keyword_set(newdarks) then begin
         gd = where(rawdarks ne 'identifyi_raw.fits')
         rawdarks=rawdarks[gd]
      endif
   endfor

   for idy=0, n_elements(rawdarks)-1 do begin
      ;; We just want the actual fits files
      ;rawdarks[idy] = strmid(rawdarks[idy],10,9) + '_raw.fits'
      ; This would be just the original names. But we want the blv files
                                ;rawdarks[idy] = strmid(rawdarks[idy],8,24)
      
      if keyword_set(precte) then begin
                                ; check directory exists. if it does
                                ; not, make the directory and copy files
         prectename = precte_dataloc+darkfiles_root[idx]+'_drk'
         checka = FILE_TEST(prectename)
         if checka eq 0 then FILE_MKDIR, prectename
         checkfile = FILE_TEST(precte_dataloc + darkfiles_root[idx]+'_drk/'+rawdarks[idy])
         if checkfile eq 0 then spawn, 'cp -a ' + rawdarks[idy] + ' ' + precte_dataloc + darkfiles_root[idx]+'_drk/'
      endif else begin
         postctename = orig_dataloc+darkfiles_root[idx]+'_drk'
         checkb = FILE_TEST(postctename)
         if checkb eq 0 then FILE_MKDIR, postctename
         checkfile = FILE_TEST(orig_dataloc + darkfiles_root[idx]+'_drk/'+rawdarks[idy])
         if checkfile eq 0 then spawn, 'cp -a ' + rawdarks[idy] + ' ' + orig_dataloc + darkfiles_root[idx]+'_drk/'

      endelse

;else begin
;      postctename = orig_dataloc+darkfiles_root+'_drk'
;      checkb = FILE_TEST(postctename)
;      if checkb eq 0 then FILE_MKDIR, postctename
;      checkfile = FILE_TEST(orig_dataloc + darkfiles_root+'_drk/'+rawdarks[idy])
;      if checkfile eq 0 then spawn, 'cp -a ' + rawdarks[idy] + ' ' + orig_dataloc + darkfiles_root+'_drk/'
;   endelse
   
   endfor
   

   ;if keyword_set(ctecor) then begin
   ;   prectename = precte_dataloc+darkfiles_root[idx]+'_drk'
   ;   checka = FILE_TEST(prectename)
   ;   ; Note: here we check to make sure it exists before we gzip
   ;   if checka eq 1 then spawn, 'gzip ' + precte_dataloc + darkfiles_root[idx]+'_drk/*.fits'
   ;   
   ;   ; Ok - now lets run cte correction code
   ;   cd, orig_dataloc + darkfiles_root[idx]+'_drk/'
   ;   runctecorr, /mp
   ;   ractoraw
   ;   cd, inputdatadir
   ;endif

endfor

cd, origdir


end
