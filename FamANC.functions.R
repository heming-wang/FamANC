#### Source functions for FamANC analyses
## Genetic Distance file
Genetic.Distance=function(map,gmap){
  cm=as.numeric()
  for(j in 1:dim(map)[1]){
    dis=gmap[,2]-map[j,4]
    point=which(dis>=0)[1]
    if(dis[point]==0){cm[j]=gmap[point,4]}
    ### 1st column bp, 4th column cm
    else{
      x1=gmap[(point-1),2]
      x2=gmap[point,2]
      x=map[j,4]
      y1=gmap[(point-1),4]
      y2=gmap[point,4]
      cm[j]=(y2-y1)*(x-x1)/(x2-x1)+y1
    }
  }

  Morgan=cm/100
  return (Morgan)


}


## Indentify possible Mendelian paths given fam data.
library("kinship2")

Mendelian.Path=function(ped){

  names(ped)=c("FID","IID","PAT","MAT","SEX","PHE")
  ped$no=1:dim(ped)[1]

  ped$depth=kindepth(id=ped$IID,dad.id=ped$PAT,mom.id=ped$MAT)

  for(i in 1:dim(ped)[1]){
    ped$parent.no[i] = sum(c(ped$PAT[i],ped$MAT[i]) %in% ped$IID)
  }
  ped=ped[order(10*ped$depth+ped$parent.no),]

  no.founder=sum(ped$parent.no==0)
  MM=expand.grid(rep(list(0:2),no.founder))

  tr2=matrix(c(0.5,0.5,0,
               0.25,0.5,0.25,
               0,0.5,0.5),nrow=3,byrow=T)
  dimnames(tr2)=list(parent=c(0,1,2),offspring=c(0,1,2))

  tr3=matrix(c(1,0,0,
               0.5,0.5,0,
               0,1,0,
               0.5,0.5,0,
               0.25,0.5,0.25,
               0,0.5,0.5,
               0,1,0,
               0,0.5,0.5,
               0,0,1),nrow=9,byrow=T)
  dimnames(tr3)=list(parents=c("00","01","02","10","11","12","20","21","22"),offspring=c(0,1,2))


  for(i in (no.founder+1):dim(ped)[1]){
    m1=rep(0:2,each=dim(MM)[1])
    MM=rbind(MM,MM,MM)
    MM=cbind(MM,m1)
    if(ped$parent.no[i]==1){
      parent=which(ped$IID %in% c(ped$PAT[i],ped$MAT[i]))
      MM=MM[-which(abs(MM[,parent]-MM[,i])==2),]
    } else {
      father=which(ped$IID == ped$PAT[i])
      mother=which(ped$IID == ped$MAT[i])
      p=apply(MM,1, function(g)  tr3[paste(g[father],g[mother],sep=""),paste(g[i])])
      MM=MM[p!=0,]
    }

  }

  names(MM)=ped$IID
  return(list(ped.path=MM,ped=ped))
}


## Estimate error probability between observed path and Mendelian path

Prob.Error=function(g1,g2,epsilon){
  Pe=matrix(c((1-epsilon)^2,2*epsilon*(1-epsilon),epsilon^2,
              epsilon*(1-epsilon),(1-epsilon)^2+epsilon^2,epsilon*(1-epsilon),
              epsilon^2,2*epsilon*(1-epsilon),(1-epsilon)^2),nrow=3,byrow=T)
  dimnames(Pe)=list(true=c(0,1,2),obs=c(0,1,2))
  p=1
  for(i in 1:length(g1)){
    p=p*Pe[paste(g1[i]),paste(g2[i])]
  }
  return(p)
}


## Choose optimal path
Optimal.Path=function(ped,g,ped.path,thres,epsilon){


  err=apply(ped.path,1,function(x) sum(x!=g))
  if(min(err==0)){
    optimal.path=g
    pr=1
  } else {
    #min.err.path=ped.path$path[err==min(err),]
    #pr=ped.path$probablity[err==min(err)]
    min.err.path=ped.path[err> 0 & err ,]
    pr=apply(min.err.path,1,function(x) Prob.Error(x,g,epsilon))
    #pr=ped.path$probablity[err<=thres]
    #pr=pr/sum(pr)

  }

  return(list(optimal.path=min.err.path,pr=pr))
}






##
Pedigree.Hap=function(anc){
  hap=apply(anc,2,function(x) paste(x,sep="",collapse =""))
  z=rle(hap)
  hap=sapply(z$values,function(x) as.numeric(strsplit(x,"")[[1]]))
  return(list(length=z$lengths,haplotype=hap))
}


## Correcting ancestry

Mendelian.Anc=function(anc,mPath,Morgan,t=6,thres=2,epsilon){

  ped.path=mPath$ped.path
  anc=anc[mPath$ped$no,]
  ped=mPath$ped

  pie=rowMeans(anc)/2
  hap=apply(anc,2,function(x) paste(x,sep="",collapse =""))
  z=rle(hap)

  Morgan.l=Morgan.r=as.numeric()
  for(i in 2:(length(z$values)-1)){
    sep=sum(z$lengths[1:i])
    Morgan.r[i-1]=Morgan[sep]
    Morgan.l[i-1]=Morgan[sep+1]
  }

  pool=apply(ped.path,1,function(x) paste(x,sep="",collapse =""))


  logic1=sapply(z$values,function(x) x %in% pool)
  logic2=z$lengths>20
  logic=logic1 & logic2


  haplotype.mendel=haplotype=sapply(z$values,function(x) as.numeric(strsplit(x,"")[[1]]))

  for(i in 1:dim(haplotype)[2]){
    if(!logic[i]){

      path = Optimal.Path(ped,haplotype[,i],ped.path,thres,epsilon)
      if(length(path$pr)==1){
        haplotype.mendel[,i]=haplotype[,i]
      } else {
        n1=max(which(logic[1:(i-1)]))
        n2=min(which(logic[(i+1):dim(haplotype)[2]]))+i
        H1=haplotype[,n1]
        H2=haplotype[,n2]
        H=path$optimal.path
        pH=rep(1,dim(H)[1])
        d1=Morgan.l[i-1]-Morgan.r[n1]
        d2=Morgan.l[n2-1]-Morgan.r[i]

        for(j in 1:dim(H)[1]){
          for(k in 1:dim(H)[2]){
            for(d in 1:2){
              p00=exp(-d*t)+(1-pie[k])*(1-exp(-d*t))
              p01=pie[k]*(1-exp(-d*t))
              p10=(1-pie[k])*(1-exp(-d*t))
              p11=exp(-d*t)+pie[k]*(1-exp(-d*t))
              transition.prob=matrix(c(p00*p00,p01*p00+p00*p01,p01*p01,
                                       p10*p00,p10*p01+p00*p11,p01*p11,## divided by 2
                                       p10*p10,p11*p10+p10*p11,p11*p11),nrow=3,byrow=T)

              dimnames(transition.prob)=list(state=c(0,1,2),state=c(0,1,2))
              pH[j]=pH[j]*transition.prob[paste(H1[k]),paste(H[j,k])]*transition.prob[paste(H[j,k]),paste(H2[k])]
            }





          }
          pH[j]=path$pr[j]*pH[j]
        }
        haplotype.mendel[,i]=t(H[which.max(pH),])

      }

    }

  }


  anc.mendel=t(apply(haplotype.mendel,1,function(x) rep(x,times=z$lengths)))

  return(anc.mendel)

}


## plot pedigree
plotped=function(ped){
  names(ped)=c("FID","IID","PAT","MAT","SEX","PHE")
  ped$no=1:dim(ped)[1]
  ped1=pedigree(id=ped$IID,dadid=ped$PAT,momid=ped$MAT,sex=ped$SEX,famid=ped$FID)
  plot(ped1[paste(ped$FID[1])])
}

                     
## Divide large pedigree missing first-generation into smaller pedigrees 
Pedigree.Divide=function(ped){
  names(ped)=c("FID","IID","PAT","MAT","SEX","PHE")
  ped$no=1:dim(ped)[1]
  ped$depth=kindepth(id=ped$IID,dad.id=ped$PAT,mom.id=ped$MAT)

  for(i in 1:dim(ped)[1]){
    ped$parent.no[i] = sum(c(ped$PAT[i],ped$MAT[i]) %in% ped$IID)
  }
  
  ped=ped[order(10*ped$depth+ped$parent.no),]
  ped$PAT[ped$depth==0]=0
  ped$MAT[ped$depth==0]=0
  
  FID1=c(1:dim(ped)[1])

  no.founder=sum(ped$parent.no==0)
  if(no.founder<dim(ped)[1]){
    group=list()
    for(i in 1:dim(ped)[1]){
      group[[i]]=i
    }
    
    for(i in (no.founder+1):dim(ped)[1]){
      
      if(ped$parent.no[i]==1){
        parent=which(ped$IID %in% c(ped$PAT[i],ped$MAT[i]))
        FamID=FID1[c(parent,i)]
        group[[parent]]=group[[i]]=c(group[[parent]],group[[i]])
        FID1[group[[i]]]=FID1[parent]
        
      } else {
        father=which(ped$IID == ped$PAT[i])
        mother=which(ped$IID == ped$MAT[i])
        FamID=FID1[c(father,mother,i)]
        group[[father]]=group[[mother]]=group[[i]]=c(group[[father]],group[[mother]],group[[i]])
        
        
        FID1[group[[i]]]=min(FID1[c(father,mother)])
        
        
      }
    }
    for(i in (no.founder+1):dim(ped)[1]){
      
      if(ped$parent.no[i]==1){
        parent=which(ped$IID %in% c(ped$PAT[i],ped$MAT[i]))
        FamID=FID1[c(parent,i)]
        group[[parent]]=group[[i]]=c(group[[parent]],group[[i]])
        FID1[group[[i]]]=FID1[parent]
        
      } else {
        father=which(ped$IID == ped$PAT[i])
        mother=which(ped$IID == ped$MAT[i])
        FamID=FID1[c(father,mother,i)]
        group[[father]]=group[[mother]]=group[[i]]=c(group[[father]],group[[mother]],group[[i]])
        
        
        FID1[group[[i]]]=min(FID1[c(father,mother)])
        
        
      }
    }
    
    
  }


  ped$FID=paste(ped$FID,"_",FID1,sep="")
  ped=ped[order(ped$FID),]
  return(ped[,c(1:7)])

}

