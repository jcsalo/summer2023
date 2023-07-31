


table_nom<-function(Vars,vars_desc,Subset=NA, Group="Race",df=la_collect){
  df_f<-df
  if(Subset=='Adeno'){
    df_f<-df[df$Histo=="Adeno",]
    Cohort<-"Adenocarcinoma Cohort"
  } else if (Subset=='SCCA'){
    df_f<-df[df$Histo=="SCCA",]
    Cohort<-"Squamous Cell Cohort"
  } else {
    df_f<-df
    Cohort<-"Whole Cohort"
  }
  Vars_f<-df_f[,Vars]
  Group_f<-df_f[,Group]
  
  Cap = paste(Cohort,vars_desc,"grouped by",Group)
  
  tableNominal(vars=Vars_f, group=Group_f,cumsum=FALSE,cap=Cap,print.pval="fisher")
  
}




table_cont<-function(Vars,vars_desc,Subset=NA, Group="Race",df=la_collect){
  df_f<-df
  if(Subset=='Adeno'){
    df_f<-df[df$Histo=="Adeno",]
    Cohort<-"Adenocarcinoma Cohort"
  } else if (Subset=='SCCA'){
    df_f<-df[df$Histo=="SCCA",]
    Cohort<-"Squamous Cell Cohort"
  } else {
    df_f<-df
    Cohort<-"Whole Cohort"
  }
  Vars_f<-df_f[,Vars]
  Group_f<-df_f[,Group]
  
  Cap = paste(Cohort,vars_desc,"grouped by",Group)
  
  tableContinuous(vars=Vars_f, group=Group_f,stats=c("n","median"),cap=Cap,print.pval="kruskal")
}




function.uniglms<-function(outcome,varslist,glm_family,glm_dataframe)
{
  uni_fit<-""
  for(vars in varslist){
    uni_fit1<-glm(outcome ~get(vars), family=glm_family, data=glm_dataframe)
    uni_fit_tidy<-tidy(uni_fit1)
    uni_fit_tidy<-uni_fit_tidy[-1,]
    # Added Marhc 5, 2020
    uni_fit_tidy$odds_ratio<-exp(uni_fit_tidy$estimate)
    uni_fit_tidy$CI25<-exp(uni_fit_tidy$estimate + qnorm(0.025) * uni_fit_tidy$std.error)
    uni_fit_tidy$CI975<-exp(uni_fit_tidy$estimate + qnorm(0.975) * uni_fit_tidy$std.error)
    
    
    
    uni_fit_tidy[c(2:8)]<-lapply(uni_fit_tidy[c(2:8)], function(x) signif(x,digits=5) )
    col1<-paste(vars,substr(uni_fit_tidy[,1],10,100),sep="|")
    #col1<-vars
    uni_fit_tidy[,1]<-col1
    uni_fit<-rbind(uni_fit,uni_fit_tidy)
  }
  uni_fit$sig<-""
  uni_fit$sig[as.numeric(as.character(uni_fit$p.value))<0.15]<-"*"
  row.names(uni_fit)<-seq_along(1:dim(uni_fit)[1])
  uni_fit$sig[uni_fit$term=='']<-"--"
  uni_fit<-uni_fit[,c(1,5,6,7,8,9)]
  return(uni_fit)
}

function.multiglm<-function(outcome,outvar,varlist,glm_family,glm_dataframe)
{
  
  variableslist<-paste(varlist,collapse=' + ')
  form<-paste(outvar,variableslist,sep=" ~ ")
  #cat(form,'\n')
  multi_fit<-do.call("glm",list(formula=form, data=glm_dataframe, family=glm_family))
  #print(multi_fit$anova)
  #print(anova(multi_fit))
  multi_fit_STEP<-do.call("stepAIC",list(object=multi_fit,direction="both",trace=FALSE))
  #print(multi_fit_STEP$anova)
  multi_formula<-multi_fit_STEP$formula
  multi_pars<-do.call("glm",list(formula=multi_formula, data=glm_dataframe, family=glm_family))
  
  multi_pars_tidy<-tidy(multi_pars)
  multi_pars_tidy<-multi_pars_tidy[-1,]
  
  # Added Marhc 5, 2020
  multi_pars_tidy$odds_ratio<-exp(multi_pars_tidy$estimate)
  multi_pars_tidy$CI25<-exp(multi_pars_tidy$estimate + qnorm(0.025) * multi_pars_tidy$std.error)
  multi_pars_tidy$CI975<-exp(multi_pars_tidy$estimate + qnorm(0.975) * multi_pars_tidy$std.error)
  
  
  
  multi_pars_tidy[c(2:5)]<-lapply(multi_pars_tidy[c(2:5)], function(x) signif(x,digits=5) )
  multi_pars_tidy$p.value<-lapply(multi_pars_tidy$p.value, function(x) format(x,scientific=FALSE,digits=7))
  #cat("\n\nResults of final model:\n")
  #print(multi_pars_tidy, row.numbers=FALSE)
  row.names(multi_pars_tidy)<-seq_along(1:dim(multi_pars_tidy)[1])
  multi_pars_tidy<-multi_pars_tidy[,c(1,5,6,7,8)]
  return_list<-list(form,multi_pars_tidy)
  return(return_list)
}





