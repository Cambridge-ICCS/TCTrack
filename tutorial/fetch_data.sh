# Create a directory to store the data in if it does not already exist
mkdir -p data/

# Fetch data from various sources using wget and place in directory

# Most data can be obtained from ESGF CEDA portal:
# psl
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/day/psl/gn/files/d20180730/psl_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc
# sfc wind
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/day/sfcWind/gn/files/d20180730/sfcWind_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc
# orography
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/fx/orog/gn/files/d20200910/orog_fx_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn.nc
# t
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/day/ta/gn/files/d20180730/ta_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc

# u
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/day/ua/gn/files/d20180730/ua_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc
# us
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/day/uas/gn/files/d20180730/uas_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500101-19501230.nc

# v
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/day/va/gn/files/d20180730/va_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc
#vs
wget --directory-prefix data https://esgf.ceda.ac.uk/thredds/fileServer/esg_cmip6/CMIP6/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/day/vas/gn/files/d20180730/vas_day_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_19500701-19501230.nc

# zg has to be fetched direct from CEDA as not present on ESGF
wget https://dap.ceda.ac.uk/badc/cmip6/data/PRIMAVERA/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/Prim3hrPt/zg7h/gn/files/d20180730/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195008010000-195008302100.nc?download=1 -O data/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195008010000-195008302100.nc
wget https://dap.ceda.ac.uk/badc/cmip6/data/PRIMAVERA/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/Prim3hrPt/zg7h/gn/files/d20180730/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195009010000-195009302100.nc?download=1 -O data/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195010090000-195009302100.nc
wget https://dap.ceda.ac.uk/badc/cmip6/data/PRIMAVERA/HighResMIP/MOHC/HadGEM3-GC31-HM/hist-1950/r1i1p1f1/Prim3hrPt/zg7h/gn/files/d20180730/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195010010000-195010302100.nc?download=1 -O data/zg7h_Prim3hrPt_HadGEM3-GC31-HM_hist-1950_r1i1p1f1_gn_195010010000-195010302100.nc
