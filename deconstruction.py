import ddosa
import dataanalysis as da
import os
import pifs
import pyfits
from numpy import *

class CatForSpectraFromImaging(ddosa.CatForSpectraFromImaging):
    version="negsig"
    minsig=-100

class FitShadow(ddosa.ii_spectra_extract):
 #   input_gb=ddosa.ghost_bustersImage
 #   input_maps=ddosa.BinMapsImage
    input_cat=CatForSpectraFromImaging
    shdtype="BIN_S"
    binary="ii_spectra_extract"
    version="stdiis"
    #binary=os.environ["COMMON_INTEGRAL_SOFTDIR"]+"/spectral/ii_spectra_extract/ii_spectra_extract_saveflux/ii_spectra_extract"

class SpectraBins(ddosa.SpectraBins):
    run_for_hashe=True

    def main(self):
        return ddosa.ImageBins

class ii_spectra_extract(ddosa.ii_spectra_extract):
    binary="ii_spectra_extract"

    #binary=os.environ["COMMON_INTEGRAL_SOFTDIR"]+"/spectral/ii_spectra_extract/ii_spectra_extract_saveflux/ii_spectra_extract"

class Dict:
    data=None
    keys=None

    def __setattr__(self,k,v):
        if self.data is None:
            self.data

    def next(self):
        if self.data is None:
            raise Exception("empty!")

        self.data.append([])
        
    
class CleanShadow(ddosa.DataAnalysis):
    #input_gb=ddosa.ghost_bustersImage
    input_gb=ddosa.ghost_bustersSpectra
    input_maps=ddosa.BinMapsImage
    input_pif=pifs.ii_pif_fromimaging
    input_fit=FitShadow
    input_imaging=ddosa.ii_skyimage

    cached=True


    mode="BIN_S"

    version="v13"

    clean_fit_kind="fit_cppix"
    diffuse=('constant')

    pfit=True
    regress=True

    def get_version(self):
        v=self.get_signature()+"."+self.version+".cleanwith_"+self.clean_fit_kind
        for d in self.diffuse:
            v+=".diffuse"+d
        v+="pfit" if self.pfit else "nopfit"
        v+="regress" if self.regress else "noregress"
        return v

    def main(self):
        corshad_f=pyfits.open(self.input_gb.corshad.get_path())
        bkgmap_f=pyfits.open(self.input_maps.back.get_path())
        offcormap_f=pyfits.open(self.input_maps.corr.get_path())
        pifs_f=pyfits.open(self.input_pif.pifs.get_path())
        fit_f=pyfits.open(self.input_fit.spectrum.get_path())
        skyres_f=pyfits.open(self.input_imaging.skyres.get_path())

        pif_name_band=[[f.header['NAME'],f.header['E_MIN'],f.header['E_MAX']] for f in pifs_f[2:]]

        source_table=[]       
        
        n_eband=corshad_f[1].data.shape[0]/3
        for i_eband in range(n_eband):
            i=2+i_eband*3
            emin,emax=corshad_f[i].header['E_MIN'],corshad_f[i].header['E_MAX']
            print "working in energy band",i_eband,emin,emax

            corshad=corshad_f[i].data
            bkgmap=bkgmap_f[i_eband+2]
            offcormap=offcormap_f[i_eband+2]
                
            saclaynorm=128.0*128.0/2.0
            ontime=corshad_f[2+i_eband*3].header['ONTIME']
            tstart=corshad_f[2+i_eband*3].header['TSTART']

            
            n_source=len(fit_f)-2

            source_models=[]
            cfactor_nonom=1./saclaynorm*ontime

            for i_source,source_e in enumerate(fit_f[2:]):
                print
                print "source",source_e.header['NAME']

                r,e,fit_q=source_e.data['RATE'][i_eband],source_e.data['STAT_ERR'][i_eband],source_e.data['QUALITY'][i_eband]

                y,z=source_e.header['Y_FIN'],source_e.header['Z_FIN']
                iz=-y+200.
                iy=z+200.
                cnom=offcormap.data[iy,iz]
                cfactor=cfactor_nonom*cnom

                    
                print "source position used in extraction y,z; nomex ",y,z, cnom,cfactor,ontime,saclaynorm*ontime
                
                if source_e.header['NAME']=="Background":
                    cfactor/=cnom
                
                if fit_q==0 or source_e.header['NAME']=="Background":
                    fit_cppix=r*cfactor
                else:
                    print "bad fit quality, wrong units"
                    fit_cppix=r
                    r=fit_cppix/cfactor
                
                
                print "ii_spectra_extract results",r,e,fit_cppix

                # get skyres
                skyres_row=[row for row in skyres_f[i_eband+2].data if row['NAME']==source_e.header['NAME']]

                if skyres_row==[]:
                    print "no skyres"
                    skyres_r=r
                    skyres_e=e
                else:
                    skyres_r=skyres_row[0]['FLUX']
                    skyres_e=skyres_row[0]['FLUX_ERR']

                cppix=skyres_r*cfactor
                cppixe=skyres_e*cfactor
                print "skyres results",skyres_r,skyres_e,cppix,cppixe
                
                # get pif
                pif_i=[i for i,(n,e1,e2) in enumerate(pif_name_band) if n==source_e.header['NAME'] and e1==emin and e2==emax]
                if pif_i==[]:
                    print "no pif found!"
                    if source_e.header['NAME']=="Background":
                        pif_e=bkgmap
                        pif_e.data/=average(pif_e.data)
                        print "will use background",bkgmap.data.sum()
                    else:
                        continue
                else:
                    pif_e=pifs_f[pif_i[0]+2]
                    assert(pif_e.header['NAME']==source_e.header['NAME'])


                   #corshad_f[2+i_eband*3].data-=pif_fitted*corshad_f[2+i_eband*3+2].data

                pif_fitted=pif_e.data*fit_cppix
                print "  fitted pif sum",pif_fitted.sum()
                #pif_fitted[isnan(pif_fitted)]=0 #??
                #sum_source_model+=pif_fitted
                
                source_models.append({'source_name':source_e.header['NAME'],'pif':pif_e.data,'fit_cppix':fit_cppix,'cfactor':cfactor,'fit_sigma':r/e,'skyres_cppix':cppix, \
                                      'fit_r':r,'fit_e':e,'skyres_r':skyres_r,'skyres_e':skyres_e,'nomex':cnom,'y':y,'z':z})
                # memory?

            # add more models
                
            if "constant" in self.diffuse:
                source_models.append({'source_name':'constant','pif':pif_e.data*0+1.,'fit_cppix':1,'cfactor':cfactor_nonom,'fit_sigma':1,'skyres_r':0,'skyres_e':1,'fit_r':1,'fit_e':1,'y':0,'z':0,'nomex':0,'skyres_cppix':0})

            def plane_model(dy,dz):
                y,z=meshgrid(arange(130),arange(134))
                return (y-65)*dy+(z-67)*dz
                
            if "plane" in self.diffuse:
                source_models.append({'source_name':'plane','pif':((0,0),plane_model),'fit_cppix':1,'cfactor':cfactor_nonom,'fit_sigma':1,'skyres_r':0,'skyres_e':1,'fit_r':1,'fit_e':1,'y':0,'z':0,'nomex':0,'skyres_cppix':0})

            #//

            if self.pfit:
                self.fit_pif(source_models,corshad_f[2+i_eband*3].data,corshad_f[2+i_eband*3+1].data,corshad_f[2+i_eband*3+2].data)
            else:
                for s in source_models:
                    s['pfit_cppix']=0
            
            if self.regress:
                self.regress_pif(source_models,corshad_f[2+i_eband*3].data,corshad_f[2+i_eband*3+1].data,corshad_f[2+i_eband*3+2].data)
            else:
                for s in source_models:
                    s['regress_cppix']=0
                    s['regress_cppix_err']=0

            for kind in ['fit_cppix','skyres_cppix']+ (['regress_cppix'] if self.regress else []) + (['pfit_cppix'] if self.pfit else []):
                self.evaluate_residuals(kind,source_models,corshad_f[2+i_eband*3].data,corshad_f[2+i_eband*3+1].data,corshad_f[2+i_eband*3+2].data)
                
            # save
            for a in source_models:
                source_table.append([a['source_name'],tstart,emin,emax,ontime,a['skyres_r'],a['skyres_e'],a['fit_r'],a['fit_e'],a['y'],a['z'],a['nomex'],a['skyres_cppix'],a['fit_cppix'],a['pfit_cppix'],a['regress_cppix'],a['regress_cppix_err']])

            # combine_model
        
            sum_source_model=zeros_like(corshad_f[2].data)  
            for a in source_models:
                if isinstance(a['pif'],tuple):
                    sum_source_model+=a[self.clean_fit_kind]*a['pfit_model']
                else:
                    sum_source_model+=a[self.clean_fit_kind]*a['pif']
            
            print "model",sum_source_model.min(),sum_source_model.max()

            print "before:",corshad_f[2+i_eband*3].data.min(),corshad_f[2+i_eband*3].data.max()
            if self.mode=="BIN_S":
                corshad_f[2+i_eband*3].data=corshad_f[2+i_eband*3].data/corshad_f[2+i_eband*3+2].data-sum_source_model
                corshad_f[2+i_eband*3].data[isnan(corshad_f[2+i_eband*3].data) | isinf(corshad_f[2+i_eband*3].data)]=0 #?
                corshad_f[2+i_eband*3].header['ISDCLEVL']="BIN_I"
                corshad_f[1].data['ISDCLEVL']="BIN_I"
            elif self.mode=="BIN_I":
                corshad_f[2+i_eband*3].data-=sum_source_model
            else:
                raise
            print "after:",corshad_f[2+i_eband*3].data.min(),corshad_f[2+i_eband*3].data.max()
            corshad_f[2+i_eband*3].data-=corshad_f[2+i_eband*3].data.min()


        savetxt("source_results.txt",array([[str(e).replace(" ","_") for e in l] for l in source_table]),fmt="%s")
        self.source_results=da.DataFile("source_results.txt")

        fn="isgri_cor_shad_cleaned.fits"
        corshad_f.writeto(fn,clobber=True)

        self.corshad=da.DataFile(fn)
                
 #           pifs=offcormap_f[i_eband+2]

    #def get_pif():
    
    def evaluate_residuals(self,kind,source_models,data,data_variance,efficiency):
        residuals=copy(data)
        for a in source_models:
            if isinstance(a['pif'],tuple):
                print "can not handle this.."
            else:
                residuals-=a[kind]*a['pif']*efficiency
        r=(residuals**2/data_variance)
        m=~(isnan(r) | isinf(r))
        rs=r[m].sum()
        ndof=where(m)[0].shape[0]
        print "residuals for",kind,rs,rs/ndof,ndof
    
    def regress_pif(self,source_models,data,data_variance,efficiency):
        import statsmodels.api as sm

        models=[]
        m=efficiency>0.1
        for a in source_models:
            if isinstance(a['pif'],tuple):
                print "can not handle this.."
            else:
                models.append((a['pif']*efficiency)[m])

        mod_wls = sm.WLS(data[m], transpose(models), weights=data_variance[m]**-0.5)
        res_wls = mod_wls.fit()
        p_best=res_wls.params
        p_best_err=res_wls.bse
        
        print(res_wls.summary())

        i=0
        for a in source_models:
            if isinstance(a['pif'],tuple):
                print "can not handle this.."
                a['regress_cppix']=0
                a['regress_cppix_err']=0
            else:
                print a['source_name'],a['fit_cppix'],a['fit_cppix']/a['fit_sigma'],p_best[i],p_best_err[i]
                a['regress_cppix']=p_best[i]
                a['regress_cppix_err']=p_best_err[i]
                i+=1

    def regress_pif_sklearn(self,source_models,data,data_variance,efficiency):
        models=[]
        for a in source_models:
            if isinstance(a['pif'],tuple):
                print "can not handle this.."
            else:
                models.append((a['pif']*efficiency).flatten())

        from sklearn import linear_model
        clf = linear_model.Ridge(fit_intercept=True, normalize=False)
        clf.fit(transpose(models),data.flatten(),sample_weights=data_variance**-0.5)
        p_best=clf.coef_

        
        i=0
        for a in source_models:
            if isinstance(a['pif'],tuple):
                print "can not handle this.."
            else:
                print a['source_name'],a['fit_cppix'],a['fit_cppix']/a['fit_sigma'],p_best[i]
                a['regress_cppix']=p_best[i]
                i+=1

    def fit_pif(self,source_models,data,data_variance,efficiency):

        p0=[]
        for a in source_models:
            p0.append(a['fit_cppix'])
            if isinstance(a['pif'],tuple):
                p0+=list(a['pif'][0])

        def residual_func(p,gp):
            s=None
            i=0
            for a in source_models:
                if isinstance(a['pif'],tuple):
                    model=p[i]*a['pif'][1](*p[i+1:i+1+len(a['pif'][0])])
                    i+=1+len(a['pif'][0])
                else:
                    model=p[i]*a['pif'] # a duck
                    i+=1
                if s is None: s=zeros_like(model)
                s+=model

            r=((s*efficiency-data)**2/data_variance)
            rs=r[(efficiency>0.1) & (data_variance>0)].sum()
            if isnan(rs): return float('inf')
            return rs.astype(float64)
    
        import nlopt
        print "p0",p0

        #opt=nlopt.opt(nlopt.LN_PRAXIS, len(p0))
        #opt=nlopt.opt(nlopt.LN_COBYLA, len(p0))
        opt=nlopt.opt(nlopt.LN_BOBYQA, len(p0))
        opt.set_lower_bounds([-1000]*len(p0))
        opt.set_upper_bounds([1000]*len(p0))
        opt.set_min_objective(residual_func)
        opt.set_xtol_rel(1e-6)
        p_best=opt.optimize(p0)
        
        i=0
        for a in source_models:
            print a['source_name'],a['fit_cppix'],a['fit_cppix']/a['fit_sigma'],p_best[i]
            a['pfit_cppix']=p_best[i]
            if isinstance(a['pif'],tuple):
                model=p_best[i]*a['pif'][1](*p_best[i+1:i+1+len(a['pif'][0])])
                a['pfit_model']=model
                print "...",p_best[i+1:i+1+len(a['pif'][0])],a['pif'][0],i,len(p_best)
                pyfits.PrimaryHDU(model).writeto(a['source_name']+".fits",clobber=True)
                fn=a['source_name']+".txt"
                savetxt(fn,p_best[i:i+1+len(a['pif'][0])])
                setattr(self,fn,da.DataFile(fn))
                i+=1+len(a['pif'][0])
            else:
                i+=1

class VerifyMDUEfficiency(ddosa.DataAnalysis):
    input_gb=ddosa.ghost_bustersImage
    input_cleaned=CleanShadow
    mode="BIN_I" # !!

    cached=True

    version="v2"

    copy_cached_input=False

    def main(self):
        corshad_f=pyfits.open(self.input_gb.corshad.get_path())
        corshad_cleaned_f=pyfits.open(self.input_cleaned.corshad.get_path())
            
        i,j=meshgrid(arange(128),arange(128))

        ei,ej=i+(i/64)*2,j+(j/32)*2
        mdu=((127-j)/32)+((127-i)/64)*4
        emdu=zeros((134,130))-1
        emdu[ej,ei]=mdu[j,i]

        pyfits.PrimaryHDU(mdu).writeto("mdu.fits",clobber=True)

        mduresults=[]
        
        n_eband=corshad_f[1].data.shape[0]/3
        for i_eband in range(n_eband)[:1]:
            i_e=2+i_eband*3
    
            emin,emax=corshad_f[i_e].header['E_MIN'],corshad_f[i_e].header['E_MAX']
            print i_eband,emin,emax
                
            corshad=corshad_f[i_e].data
            corshad_var=corshad_f[i_e+1].data
            corshad_effi=corshad_f[i_e+2].data
            corshad_cleaned=corshad_cleaned_f[i_e].data

            #corshad-=average(corshad[corshad_effi>0])
           # pyfits.PrimaryHDU(corshad).writeto("corshad_nooffset.fits",clobber=True)

            mducorr=[]
            for mdu_i in range(8):
                corshad=corshad_f[i_e].data[emdu==mdu_i]
                corshad_var=corshad_f[i_e+1].data[emdu==mdu_i]
                corshad_effi=corshad_f[i_e+2].data[emdu==mdu_i]
                corshad_cleaned=corshad_cleaned_f[i_e].data[emdu==mdu_i]

                mask=corshad_effi>0.1
                mask_area=mask.shape[0]
                av=(corshad_cleaned/corshad_var)[mask].sum()/(1/corshad_var)[mask].sum()
                av_err=1/(1/corshad_var**2)[mask].sum()**0.5
                print "corshad",mdu_i,corshad[mask].sum()/mask_area,corshad_cleaned[mask].sum()/mask_area,mask_area,av,av_err
                
                t=zeros_like(corshad_f[i_e].data)
                #t[emdu==mdu_i]=corshad_cleaned-corshad_cleaned[mask].sum()/mask_area
                #t[emdu==mdu_i][mask]=-1000
                #pyfits.PrimaryHDU(t).writeto("mdu_%i.fits"%mdu_i,clobber=True)


                mduresults.append([i_eband,emin,emax,mdu_i,corshad[mask].sum()/mask_area,av,av_err])

            
            data=array(mduresults)[:,4]
            residual=array(mduresults)[:,5]
            residual_err=array(mduresults)[:,6]
            offset=average(residual)
            print "average offset",offset,(residual-offset)/residual,(residual_err)/residual

            mducorr=(residual-offset)/residual
            
            for mdu_i in range(8):
                corshad_cleaned_f[i_e].data[emdu==mdu_i]/=residual[mdu_i]/offset
                corshad_f[i_e].data[emdu==mdu_i]/=residual[mdu_i]/offset
                mduresults[mdu_i].append(residual[mdu_i]/offset)
                mduresults[mdu_i].append(residual_err[mdu_i]/offset)
                #corshad_cleaned_f[i_e].data[emdu==mdu_i]-=residual[mdu_i]
            #corshad_cleaned_f[i_e].data-=corshad_cleaned_f[i_e].data.min()

        fn="isgri_cor_shad_mduadjusted.fits"
        corshad_f.writeto(fn,clobber=True)
        self.corshad=da.DataFile(fn)
        
        fn="isgri_cor_shad_mduadjusted_cleaned.fits"
        corshad_cleaned_f.writeto(fn,clobber=True)
        self.corshad_cleaned=da.DataFile(fn)
                
        savetxt("mdu_results.txt",mduresults)
        self.mduresults=da.DataFile("mdu_results.txt")

class VerifyMDUEfficiencyP2(VerifyMDUEfficiency):
    input_gb=ddosa.ghost_bustersImage
    input_cleaned=VerifyMDUEfficiency

            
class ii_skyimage_I3(ddosa.ii_skyimage):
    input_gb=VerifyMDUEfficiency
    image_tag="I3"

class ii_skyimage(ddosa.ii_skyimage):
    image_tag="std"

class ii_skyimage_I2(ddosa.ii_skyimage):
    input_gb=CleanShadow
    image_tag="I2"

    version="v4"

    def post_process(self):
        ima=pyfits.open(self.skyima.get_path())
        neband=(len(ima)-2)/5
        table=[]
        for ieband in range(neband):
            sigh=ima[2+ieband*5+2]
            exph=ima[2+ieband*5+4]
            noise_total=std(sigh.data)
            noise_center=std(sigh.data[exph.data>average(exph.data)*0.3])
            max_center=(sigh.data[exph.data>average(exph.data)*0.3]).max()
            min_center=(sigh.data[exph.data>average(exph.data)*0.3]).min()

            table.append([sigh.header['E_MIN'],sigh.header['E_MAX'],noise_total,noise_center,max_center,min_center])
            print table[-1]

        savetxt("image_stats.txt",table)
        self.imastats=da.DataFile("image_stats.txt")

class FitShadowI2(FitShadow):
    input_gb=CleanShadow
    version="v2"
    shdtype="BIN_S"

