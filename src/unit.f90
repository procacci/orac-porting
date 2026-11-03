module unit
  implicit none
  
  INTEGER,PARAMETER :: nfile_dump=20
  INTEGER           :: kread,kprint,kdump,kdump_out,kgdata,kcnfig,knlist,kgofr &
       ,kplot,kconf,kpot,kout,kfield,kvaf,krms,kxrms,kxrms_atm,ksol1 &
       ,ksol2,kgroup,ktop,ktemplate,ktest,iurest,kspec,kgr &
       ,kplot_fragm,kplot_steer,kplot_steer_temp,kplot_temp_steer &
       ,kplot_stn0,kplot_alch0,kfragm_dist,kplot_center,kgofr_sk,kvi,kcoord_slv &
       ,kcoord_slt,ktpgprm_write,ktpgprm_read,kavg,kavg_xrms &
       ,kwrite_dump(nfile_dump),no_records(nfile_dump),nwrite_dump &
       ,ktopol,khbonds,kabmd,kdipole,kdiff,kfreq,kdynamic,kwork,kitp
  CHARACTER(LEN=80) ::  file_names(nfile_dump)
  REAL(8)           :: pi,efact,avogad,boltz,gascon,planck,elechg,epso,twopi,boxl
  REAL(8)           :: unitm,unitl,unitt,uniti,unite,unitc,unitp,ttstep
  character*10       :: ffield_type
  
end module unit
