library(readxl)
library(gwasrapidd)
library(dplyr)
library(ggpubr)
options(bitmapType='cairo')
shortenname<-function(xx,L){
    res=do.call(c,lapply(xx,function(x){
    if (nchar(x)>L) {
	    paste0(substr(x,1,40),'...')
    }else x
    }))
    return (res)
}
local_clump<-function(exp_dat,clump_p=1,clump_kb=10000,clump_r2=0.001){
        fitld<-try(exposure <-ld_clump_local(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, 
                id=exp_dat$exposure),clump_p=clump_p,clump_kb=clump_kb,clump_r2=clump_r2,
                bfile="1kg.v3/EUR", #https://github.com/MRCIEU/gwasvcf
                #bfile='loaddat//g1000_eur/EUR',
                plink_bin=get_plink_exe()),silent=FALSE)
        if ('try-error' %in% class(fitld)){
            return (exp_dat)
        }
        return (exp_dat[exp_dat$SNP %in% exposure$rsid,])
}


MR_gut<-function(outcome_id=finn-b-L12_HYPETROPHICSCAR,pvalue=1e-5,method=c('mr_ivw'),output=c('html'),outdir='./'){
	print('this is the latest version')
	method=head(method,5)
	pvalue=as.numeric(pvalue)

	final=c()
	Fstat_sig=c()
	presso_res=c()
	heter_res=c()
	pleio_res=c()
	direc_res=c()
	tmpdir=gsub('result','tmp',dirname(outdir))
	if (file.exists(paste0(tmpdir,'/mr_result.rda'))) load(paste0(tmpdir,'/mr_result.rda'))
	if (file.exists(paste0(tmpdir,'/Fstat_sig.rda'))) load(paste0(tmpdir,'/Fstat_sig.rda'))
	if (file.exists(paste0(tmpdir,'/presso_res.rda'))) load(paste0(tmpdir,'/presso_res.rda'))
	if (file.exists(paste0(tmpdir,'/heter_res.rda'))) load(paste0(tmpdir,'/heter_res.rda'))
	if (file.exists(paste0(tmpdir,'/pleio_res.rda'))) load(paste0(tmpdir,'/pleio_res.rda'))
	if (file.exists(paste0(tmpdir,'/direc_res.rda'))) load(paste0(tmpdir,'/direc_res.rda'))
	for (id in dir('gut')){
		#if (as.integer(system(paste0("grep ",id,' ',tmpdir,'/MR_',dir2load,'_done|wc -l'),intern=T))>0) next
		exposure=read.delim(paste0('gut/',id),header=T,sep='\t')
		names(exposure)=c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure','se.exposure','pval.exposure','samplesize.exposure','exposure','eaf.exposure'
,'id.exposure')
		if (file.exists('1kg.v3')) {
			exposure=local_clump(exposure,clump_kb=10000,clump_r2=0.001)
			}else exposure=TwoSampleMR::clump_data(exposure,clump_kb=10000,clump_r2=0.001)

		cat(id,' is running','\n')
		exposure=exposure[exposure$pval<pvalue,]
		if (!is.null(outfile)) outcome=read.delim(outfile,header=T,sep='\t')
		if (!is.null(outcome_id)) outcome=TwoSampleMR::extract_outcome_data(snps=unique(exposure$SNP),outcomes=outcome_id,access_token=NULL)
		outcome=outcome[outcome$SNP %in% exposure$SNP,]
		outcome=outcome[,c('SNP','effect_allele.outcome','other_allele.outcome','beta.outcome','se.outcome','pval.outcome','samplesize.outcome','outcome','eaf.outcome'
,'id.outcome')]
		if (dim(exposure)[1]==0)next
		#merging and rm incompatable SNPs
		#exposure$id.exposure=exposure$exposure
	        dat=merge(exposure,outcome,by='SNP')
        	dat=dat[dat$effect_allele.exposure==dat$effect_allele.outcome & dat$other_allele.exposure==dat$other_allele.outcome,]
	        dat=TwoSampleMR::harmonise_data(exposure_dat=exposure,outcome_dat=outcome)
        	dat=unique(dat)
	        if (dim(dat)[1]<3) next
	        dat[,'r.exposure']=TwoSampleMR::get_r_from_bsen(dat$beta.exposure,dat$se.exposure,max(exposure$samplesize.exposure))**2
        	dat[,'r.outcome']=TwoSampleMR::get_r_from_bsen(dat$beta.outcome,dat$se.outcome,max(outcome$samplesize.outcome))**2
		dat$r.exposure[is.na(dat$r.exposure)]=0
		dat$r.outcome[is.na(dat$r.outcome)]=0
	        #write.table(dat,paste0(outdir,'/MR_harmonise.txt'),row.names=F,col.names=T,quote=F,sep='\t')
        	#MR procedure
	        #proxy=paste(unlist(lapply(unique(dat$exposure),function(x) x)),collapse='; ')  
        	#dat$exposure=paste(unlist(lapply(unique(dat$exposure),function(x)strsplit(x,' \\|\\| ')[[1]][1])),collapse='; ')
	        #if (unique(dat$exposure)=='')dat$exposure=proxy
        	mr_res=TwoSampleMR::mr(dat,method_list=method)
		if (dim(mr_res)[1]==0) next
		if (all(is.na(mr_res$pval)))next
		if (is.na(mr_res$pval[mr_res$method=='Inverse variance weighted'])) next
	        mr_result=TwoSampleMR::generate_odds_ratios(mr_res)
		final=rbind(final,mr_result)
		save(final,file=paste0(tmpdir,'/mr_result.rda'))


        	if (mr_res$pval[mr_res$method=='Inverse variance weighted']<0.05) {
       			heter_res=rbind(heter_res,TwoSampleMR::mr_heterogeneity(dat,method_list=c('mr_ivw','mr_egger_regression')))
			pleiotropy=TwoSampleMR::mr_pleiotropy_test(dat)
        		pleio_res=rbind(pleio_res,pleiotropy)
       			direc_res=rbind(direc_res,TwoSampleMR::directionality_test(dat))
			suboutdir=paste0(outdir,'/',id)
			dir.create(paste0(suboutdir))
			#final=rbind(final,mr_result)
			if ('html' %in% output){
				tryCatch({
        	        		TwoSampleMR::mr_report(dat,output_path=suboutdir,output_type = "html")
	        	        },error=function(e){return (NULL)})
			}
			if ('scatter' %in% output){
				pdf(paste0(suboutdir,'/',id,'_scatter_plot.pdf'),wi=8,he=6)
			        p=TwoSampleMR::mr_scatter_plot(mr_result,dat)
			        for (i in names(p)) print(p[[i]])
			        dev.off()
			}
			if ('forest' %in% output){
				pdf(paste0(suboutdir,'/',id,'_forest_plot.pdf'),wi=6,he=8)
			        print (TwoSampleMR::mr_forest_plot(TwoSampleMR::mr_singlesnp(dat)))
			        dev.off()
			}
			if ('funnel' %in% output){
				pdf(paste0(suboutdir,'/',id,'_funnel_plot.pdf'),wi=8,he=6)
			        print(TwoSampleMR::mr_funnel_plot(TwoSampleMR::mr_singlesnp(dat)))
			        dev.off()
			}
			if ('leaveoneout' %in% output){
				pdf(paste0(suboutdir,'/',id,'_leaveOneOut_plot.pdf'),wi=6,he=8)
			        print(TwoSampleMR::mr_leaveoneout_plot(TwoSampleMR::mr_leaveoneout(dat)))
			        dev.off()
			}
			if ('MRpresso' %in% output & dim(dat)[1]>3){
				presso=MRPRESSO::mr_presso(dat,BetaOutcome='beta.outcome',BetaExposure='beta.exposure',SdOutcome='se.outcome',SdExposure='se.exposure',OUTLIERtest = TRUE, DISTORTIONtest = TRUE)
			        #presso=run_mr_presso(dat)
				presso[[1]]$id=id
				presso[[1]]$RSSobs=presso[[2]][[1]]$RSSobs
				presso[[1]]$Global.Test.Pvalue=presso[[2]][[1]]$Pvalue
			        print (presso)
				presso_res=rbind(presso_res,presso[[1]])
			        #write.table(presso,paste0(outdir,'/',id,'_presso_test.txt'),row.names=F,col.names=T,quote=F,sep='\t')

			}
			Fstat<-do.call(rbind,lapply(1:dim(exposure)[1],function(idx){
	                N=exposure$samplesize[idx]
        	        sd=exposure$se[idx]*sqrt(N)
	                #R2=2*(1-exposure$eaf.exposure[idx])*exposure$eaf.exposure[idx]*exposure$beta.exposure[idx]/sd
        	        R2=TwoSampleMR::get_r_from_bsen(exposure$beta[idx],exposure$se[idx],N)**2
                	Fvalue=(N-2)*R2/(1-R2)
	                c(sd,R2,Fvalue)
        		}))
			exposure=exposure[,names(exposure)!='data_source.exposure']
			Fstat=cbind(exposure,Fstat)
			Fstat_sig=rbind(Fstat_sig,Fstat)
			save(presso_res,file=paste0(tmpdir,'/presso_res.rda'))
			save(heter_res,file=paste0(tmpdir,'/heter_res.rda'))
			save(pleio_res,file=paste0(tmpdir,'/pleio_res.rda'))
			save(direc_res,file=paste0(tmpdir,'/direc_res.rda'))
			save(Fstat_sig,file=paste0(tmpdir,'/Fstat_sig.rda'))
		}
	        #write.table(mr_result,paste0(outdir,'/MR_result.txt'),row.names=F,col.names=T,quote=F,sep='\t',append=T)
	}
	write.table(presso_res,paste0(outdir,'/MR_presso_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(heter_res,paste0(outdir,'/heterogeneity_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(pleio_res,paste0(outdir,'/pleiotropy_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	write.table(direc_res,paste0(outdir,'/steiger_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	if (sum(final[,9]<0.05 & final$method=='Inverse variance weighted')==0){
		msg=paste0('No significant result obtained.')
                cat (msg,file=paste0(outdir,'/msg.log'))
		return ()
	}
	
	if(!is.null(Fstat_sig)){
        Fstat=data.frame(Fstat_sig)
        names(Fstat)=c(names(exposure),'sd','R2','F')
        write.table(Fstat,paste0(outdir,'/exposure_info.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	}
	mr_result=data.frame(unique(final))
	mr_result$qvalue=p.adjust(mr_result$pval,method='bonferroni')
	write.table(mr_result[order(mr_result$pval),],paste0(outdir,'/MR_result.txt'),row.names=F,col.names=T,quote=F,sep='\t')
	drawforest<-function(dt){
		dt=dt[!is.infinite(dt$or),]
		dt=subset(dt,or_uci95<100 & or_lci95>0)
		if (dim(dt)[1]<2){
			msg=paste0('Your analysis has finished.')
	                cat (msg,file=paste0(outdir,'/msg.log'))
			return ()
		}
		if (dim(dt)[1]>30)dt=dt[1:30,]
                library(forestploter)
                dt$or_uci95[dt$or_uci95>100]=dt$or[dt$or_uci95>100]*1.5
                dt$` ` <- paste(rep(" ",15),collapse=" ")
                dt$`OR (95% CI)` <- ifelse(is.na(dt$or), "",
                           sprintf("%.3f (%.3f to %.3f)",
                                   dt$or, dt$or_lci95, dt$or_uci95))

                tm<-forest_theme(base_size=12,
                        refline_col='red',
                        footnote_col='#636363',
                        footnote_fontface='italic')

                dt$pval=format(dt$pval, scientific = T,digits=2)
		dt$qvalue=format(dt$qvalue,scientific=T,digits=2)
		dt$method='mr_ivw'
                p<-forestploter::forest(
                        dt[,c(1,6,15,16,9)],
                        est=dt$or,lower=dt$or_lci95,upper=dt$or_uci95,ci_column=3,
                        xlim=c(round(min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),x_trans='log',
                        ticks_at=c(.01+round(min(dt$or_lci95)-abs(min(dt$or_lci95))*.025,2),1,round(max(dt$or_uci95)+abs(max(dt$or_uci95))*.025,2)),
                        theme=tm
                )
                pdf(paste0(outdir,'/MR_forest.pdf'),wi=10,he=max(4,dim(dt)[1]*.4))
                plot(p)
                dev.off()
                return (p)
        }
	mr_result$id.exposure=shortenname(mr_result$id.exposure,30)
	plt=drawforest(mr_result[mr_result$pval<0.05 & mr_result$method=='Inverse variance weighted',names(mr_result)!='qvalue'])
	save(plt,file=paste0(outdir,'/plt.rda'))
#	data<-reshape2::dcast(mr_result[,c('id.exposure','method','pval')],id.exposure~method)
#	CN=names(data)
#	data=data.frame(data[,-1],row.names=data[,1])
#	names(data)=CN[-1]
#	source('loaddat/MR_circle.r')
#	heatmap_circle(data,attr=NULL,showname=T,outdir=outdir)
}

#load('params_gut.rda')
#MRgut(params[['outcome_id']],params[['outcome_file']],as.numeric(params[['pvalue']]),params[['method']],params[['output']],params[['dir2load']],params[['outdir']])
