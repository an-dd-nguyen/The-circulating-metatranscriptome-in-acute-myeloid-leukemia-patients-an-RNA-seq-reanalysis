### First define a function that takes two vectors 
### (or, in the case of survival object, a matrix and
### a vector), and returns a single P-value.
### The first vector is the object we're looking for an
### association with and the second is either a logical vector
### indicating presence/absence of the bug, or numeric giving dose. We assume the order is the same
### (i.e. same elements correspond to same patient 

library(survival)
#load("Vars/age.Rdata")
#age = readRDS("../Downloads/BeatAML_statistics/Normalized/clinical/ageatdiagnosis.RDS")
#load("Vars/gender.Rdata")
#gender = readRDS("../Downloads/BeatAML_statistics/Normalized/clinical/consensus_sex.RDS")
#load("Vars/year.Rdata")
#load("Vars/diag.Rdata")

library(coin)
library(clinfun)
library(MASS)

### First NOT controlling for age/gender

getPvec<-function(var, bug){

### If too few observations, punt:
if(sum(!is.na(var))<2 | sum(!is.na(bug))<2){return(NA)}

### First handle survival
if(is.Surv(var)){

if(is.factor(bug)){
cph<-coxph(var~bug)
out<-summary(cph)[[10]][3]} else{
cph<-coxph(var~bug)
out<-summary(cph)[[7]][5]}
return(out)}

### First make it simpler by swapping var and bug if bug is categorical:

if(is.factor(bug)){
tmp<-var
var<-bug
bug<-tmp}

### Now handle factor (categorial). Two cases whether bug is logical or
### numeric

if(is.factor(var) & !is.ordered(var)){

if(is.numeric(bug)){
mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-lm(bug~var)
anv<-anova(mdl)

return(anv[[5]][1])}

if(is.logical(bug)){
mdl<-glm(bug ~ var, family = "binomial")
return(anova(mdl, test="Chisq")[2,5])}

if(is.factor(bug)){
return(chisq.test(bug, var)$p.val)}

}

if(is.ordered(var)){

if(is.numeric(bug)){
g<-match(as.character(var), levels(var))
pi<-jonckheere.test(bug, g, alt="i")$p.val
pd<-jonckheere.test(bug, g, alt="d")$p.val

if(pi<pd){return(2*pi)} else return(-pd*2)}


if(is.logical(bug)){
tst<-independence_test(var~bug)
return(sign(statistic(tst))*as.numeric(pvalue(tst)))}
}


### Now handle logical-logical

if(is.logical(var) & is.logical(bug)){
mdl<-glm(var ~ bug, family = "binomial")
return(sign(mdl[[1]][2])*anova(mdl, test="Chisq")[2,5])}

### Now handle numeric-numeric

if(is.numeric(bug) & is.numeric(var)){
mn<-min(var, na.rm=T)
if(mn<=0){var<-var+abs(mn)+1} ### Make positive
var<-log(var) ### Then log-transform
mdl<-lm(var~bug)
anv<-anova(mdl)

return(sign(mdl[[1]][2])*anv[[5]][1])}

### Finally, logical-numeric. Call the logical var and
### the numeric bug

if(is.numeric(var)){
tmp<-var
var<-bug
bug<-tmp}

mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-glm(var ~ bug, family = "binomial")
return(sign(mdl[[1]][2])*anova(mdl, test="Chisq")[2,5])

} ### end function

### Now a function to control for age
### Subset is optional vector indicating subset of samples

#getPvecA<-function(var, bug, subset=NULL)
getPvecA<-function(var, bug, age)
{
  
  ### If too few observations, punt:
  if(sum(!is.na(var))<2 | sum(!is.na(bug))<2){return(NA)}
  
  #if(!is.null(subset)){
  #age<-age[subset]
  #gender<-gender[subset]}
  
  ### Now handle survival
  if(is.Surv(var))
  {
    cph<-coxph(var~bug+age)
    o1<-summary(cph)
    cph1<-coxph(var~age)
    cph2<-coxph(var~age+bug)
    return(anova(cph1,cph2)[2,4])
  }
  
  ### First make it simpler by swapping var and bug if bug is categorical:
  
  if(is.factor(bug))
  {
    print("a")
    tmp<-var
    var<-bug
    bug<-tmp
    
  }
  
  
  ### Now handle factor (categorial). Two cases whether bug is logical or
  ### numeric
  
  if(is.ordered(var))
  {
    
    if(is.numeric(bug))
    {
      mn<-min(bug, na.rm=T)
      if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
      bug<-log(bug) ### Then log-transform
      mdl<-try(polr(var~age+bug,Hess=T))
      if(!inherits(mdl, "try-error"))
      {
        
        ctable <- coef(summary(mdl))
        pval <- pnorm(abs(ctable["bug", "t value"]), lower.tail = FALSE) * 2
        
        
        return(sign(mdl[[1]][2])*pval)
      } 
      else return(NA)
    }
    
    if(is.logical(bug) | is.factor(bug))
    {
      mdl<-try(polr(var~age+bug,Hess=T))
      if(!inherits(mdl, "try-error"))
      {
        ctable <- coef(summary(mdl))
        pval <- pnorm(abs(ctable["bugTRUE", "t value"]), lower.tail = FALSE) * 2
        
        
        return(sign(mdl[[1]][2])*pval)
      } 
      else return(NA)
    }
  }
  
  if(is.factor(var) & !is.ordered(var))
  {
    
    if(is.factor(bug)){return(NA)} # Can't really control for age/gender
    
    if(is.numeric(bug))
    {
      mn<-min(bug, na.rm=T)
      if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
      bug<-log(bug) ### Then log-transform
      mdl<-lm(bug~age+var)
      anv<-anova(mdl)
      
      return(anv[[5]][2])
    }
    
    if(is.logical(bug))
    {
      mdl<-glm(bug ~ age+var, family = "binomial")
      return(anova(mdl, test="Chisq")[3,5])
    }
  }
  
  
  ### Now handle logical-logical
  
  if(is.logical(var) & is.logical(bug))
  {
    mdl<-glm(var ~ age+bug, family = "binomial")
    return(sign(mdl[[1]][3])*anova(mdl, test="Chisq")[3,5])
  }
  
  ### Now handle numeric-numeric
  
  if(is.numeric(bug) & is.numeric(var))
  {
    #print("TRUE") 
    #print(bug)
    #print(var)
    mn<-min(var, na.rm=T)
    if(mn<=0){var<-var+abs(mn)+1} ### Make positive
    var<-log(var) ### Then log-transform
    mdl<-lm(var~ age + bug)
    anv<-anova(mdl)
    #print(mdl)
    #print(c(sign(mdl[[1]][2]), anv[[5]][1], sign(mdl[[1]][2])*anv[[5]][1]))
    return(sign(mdl[[1]]["bug"])*anv[[5]][2])
  }
  
  ### Finally, logical-numeric. Call the logical var and
  ### the numeric bug
  
  if(is.numeric(var))
  {
    #print("b")
    tmp<-var
    var<-bug
    bug<-tmp
    mn<-min(bug, na.rm=T)
    if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
    bug<-log(bug) ### Then log-transform
    mdl<-glm(var ~ age+bug, family = "binomial")
    return(sign(mdl[[1]][length(mdl[[1]])])*anova(mdl, test="Chisq")[3,5])
  }
  
  if(is.logical(var))
  {
    mn = min(bug, na.rm = T)
    if(mn <= 0){bug = bug + abs(mn) + 1}
    bug = log(bug)
    mdl = glm(var ~ age + bug, family = "binomial")
    return(sign(mdl[[1]][length(mdl[[1]])])*anova(mdl, test="Chisq")[3,5])
  }
  
  

} ### end function

### Now a function to control for age and gender and diagnosis
### Subset is optional vector indicating subset of samples

getPvecAGD<-function(var, bug, subset=NULL){

### If too few observations, punt:
if(sum(!is.na(var))<2 | sum(!is.na(bug))<2){return(NA)}

if(!is.null(subset)){
age<-age
gender<-gender
diag<-diag}

### Now handle survival
if(is.Surv(var)){
cph<-coxph(var~bug+age+gender+diag)
o1<-summary(cph)
cph1<-coxph(var~age+gender+diag)
cph2<-coxph(var~age+gender+diag+bug)
return(anova(cph1,cph2)[2,4])}

### First make it simpler by swapping var and bug if bug is categorical:

if(is.factor(bug)){
tmp<-var
var<-bug
bug<-tmp}


### Now handle factor (categorial). Two cases whether bug is logical or
### numeric

if(is.ordered(var)){

if(is.numeric(bug)){
mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-try(polr(var~age+gender+diag+bug,Hess=T))
if(!inherits(mdl, "try-error")){
ctable <- coef(summary(mdl))
pval <- pnorm(abs(ctable["bug", "t value"]), lower.tail = FALSE) * 2


return(sign(mdl[[1]]["bug"])*pval)} else return(NA)}

if(is.logical(bug) | is.factor(bug)){
mdl<-try(polr(var~age+gender+diag+bug,Hess=T))
if(!inherits(mdl, "try-error")){
ctable <- coef(summary(mdl))
pval <- pnorm(abs(ctable["bugTRUE", "t value"]), lower.tail = FALSE) * 2


return(sign(mdl[[1]]["bug"])*pval)} else return(NA)}
}

if(is.factor(var) & !is.ordered(var)){

if(is.factor(bug)){return(NA)} # Can't really control for age/gender

if(is.numeric(bug)){
mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-lm(bug~age+gender+diag+var)
anv<-anova(mdl)

return(anv[[5]][4])}

if(is.logical(bug)){
mdl<-glm(bug ~ age+gender+diag+var, family = "binomial")
return(anova(mdl, test="Chisq")["var",5])}
}


### Now handle logical-logical

if(is.logical(var) & is.logical(bug)){
mdl<-glm(var ~ age+gender+diag+bug, family = "binomial")
return(sign(mdl[[1]]["bugTRUE"])*anova(mdl, test="Chisq")[5,5])}

### Now handle numeric-numeric

if(is.numeric(bug) & is.numeric(var)){
mn<-min(var, na.rm=T)
if(mn<=0){var<-var+abs(mn)+1} ### Make positive
var<-log(var) ### Then log-transform
mdl<-lm(var~age+gender+diag+bug)
anv<-anova(mdl)

return(sign(mdl[[1]]["bug"])*anv[[5]][4])}

### Finally, logical-numeric. Call the logical var and
### the numeric bug

if(is.numeric(var)){
tmp<-var
var<-bug
bug<-tmp}

mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-glm(var ~ age+gender+diag+bug, family = "binomial")
return(sign(mdl[[1]]["bug"])*anova(mdl, test="Chisq")[5,5])

} ### end function


### Now a function to control for age and gender and year
### Subset is optional vector indicating subset of samples

getPvecAGY<-function(var, bug, age, gender, year){

### If too few observations, punt:
if(sum(!is.na(var))<2 | sum(!is.na(bug))<2){return(NA)}

#if(!is.null(subset)){
#age<-age
#gender<-gender
#year<-year}

### Now handle survival
if(is.Surv(var)){
cph<-coxph(var~bug+age+gender+year)
o1<-summary(cph)
cph1<-coxph(var~age+gender+year)
cph2<-coxph(var~age+gender+year+bug)
return(anova(cph1,cph2)[2,4])}

### First make it simpler by swapping var and bug if bug is categorical:

if(is.factor(bug)){
tmp<-var
var<-bug
bug<-tmp}


### Now handle factor (categorial). Two cases whether bug is logical or
### numeric

if(is.ordered(var)){

if(is.numeric(bug)){
mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-try(polr(var~age+gender+year+bug,Hess=T))
if(!inherits(mdl, "try-error")){
ctable <- coef(summary(mdl))
pval <- pnorm(abs(ctable["bug", "t value"]), lower.tail = FALSE) * 2


return(sign(mdl[[1]]["bug"])*pval)} else return(NA)}

if(is.logical(bug) | is.factor(bug)){
mdl<-try(polr(var~age+gender+year+bug,Hess=T))
if(!inherits(mdl, "try-error")){
ctable <- coef(summary(mdl))
pval <- pnorm(abs(ctable["bugTRUE", "t value"]), lower.tail = FALSE) * 2


return(sign(mdl[[1]]["bugTRUE"])*pval)} else return(NA)}
}

if(is.factor(var) & !is.ordered(var)){

if(is.factor(bug)){return(NA)} # Can't really control for age/gender

if(is.numeric(bug)){
mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-lm(bug~age+gender+year+var)
anv<-anova(mdl)

return(anv[[5]][4])}

if(is.logical(bug)){
mdl<-glm(bug ~ age+gender+year+var, family = "binomial")
return(anova(mdl, test="Chisq")["var",5])}
}


### Now handle logical-logical

if(is.logical(var) & is.logical(bug)){
mdl<-glm(var ~ age+gender+year+bug, family = "binomial")
return(sign(mdl[[1]]["bugTRUE"])*anova(mdl, test="Chisq")[5,5])}

### Now handle numeric-numeric

if(is.numeric(bug) & is.numeric(var)){
mn<-min(var, na.rm=T)
if(mn<=0){var<-var+abs(mn)+1} ### Make positive
var<-log(var) ### Then log-transform
mdl<-lm(var~age+gender+year+bug)
anv<-anova(mdl)

return(sign(mdl[[1]]["bug"])*anv[[5]][4])}

### Finally, logical-numeric. Call the logical var and
### the numeric bug

if(is.numeric(var)){
tmp<-var
var<-bug
bug<-tmp}

mn<-min(bug, na.rm=T)
if(mn<=0){bug<-bug+abs(mn)+1} ### Make positive
bug<-log(bug) ### Then log-transform
mdl<-glm(var ~ age+gender+year+bug, family = "binomial")
return(sign(mdl[[1]]["bug"])*anova(mdl, test="Chisq")[5,5])

} ### end function


### Now execute it

dn<-dir("Vars")
dm<-dir("Bugs")

nn<-substring(dn,first=1,last=nchar(dn)-6)
nm<-substring(dm,first=1,last=nchar(dm)-6)

Ps<-vector(length=0)
nms<-vector(length=0)
PsA<-Ps
PsAGY<-Ps
PsAGD<-Ps

for(i in 1:length(dn)){
load(paste("Vars/",dn[i],sep=""))
for(j in 1:length(dm)){
load(paste("Bugs/",dm[j],sep=""))

P<-getPvec(get(nn[i]),get(nm[j]))
PA<-getPvecA(get(nn[i]),get(nm[j]))
PAGY<-getPvecAGY(get(nn[i]),get(nm[j]))
PAGD<-getPvecAGD(get(nn[i]),get(nm[j]))
fn<-paste(nn[i],nm[j],sep=".")
Ps<-c(Ps, P)
PsA<-c(PsA, PA)
PsAGY<-c(PsAGY, PAGY)
PsAGD<-c(PsAGD, PAGD)
nms<-c(nms,fn)

##save(P, file=fn)
rm(P)
rm(PA)
rm(PAGY)
rm(PAGD)
print(fn)
#rm(get(nm[j]))
gc()}
#rm(get(nn[i]))
}

names(Ps)<-nms
save(Ps, file="Ps.Rdata")
names(PsA)<-nms
save(PsA, file="PsA.Rdata")
names(PsAGY)<-nms
save(PsAGY, file="PsAGY.Rdata")
names(PsAGD)<-nms
save(PsAGD, file="PsAGD.Rdata")

### Now do it for each diagnosis individually

load("Vars/diag.Rdata")

dn<-setdiff(dir("Vars"), "diag.Rdata")
dm<-dir("Bugs")

tbl<-table(diag)
diags<-names(tbl[tbl>20])

for(k in diags){
print(k)

w<-which(diag==k)
 
Ps<-vector(length=0)
nms<-vector(length=0)
PsA<-Ps
PsAGY<-Ps

nn<-substring(dn,first=1,last=nchar(dn)-6)
nm<-substring(dm,first=1,last=nchar(dm)-6)

for(i in 1:length(dn)){
load(paste("Vars/",dn[i],sep=""))
for(j in 1:length(dm)){
load(paste("Bugs/",dm[j],sep=""))

P<-getPvec(get(nn[i])[w],get(nm[j])[w])
PA<-getPvecA(get(nn[i])[w],get(nm[j])[w], subset=w)
PAGY<-getPvecAGY(get(nn[i])[w],get(nm[j])[w], subset=w)
fn<-paste(nn[i],nm[j],sep=".")
Ps<-c(Ps, P)
PsA<-c(PsA,PA)
PsAGY<-c(PsAGY,PAGY)
nms<-c(nms, fn)
rm(P)
rm(PA)
rm(PAGY)
print(fn)
#rm(get(nm[j]))
gc()}
#rm(get(nn[i]))
}

names(Ps)<-nms
names(PsA)<-nms
names(PsAGY)<-nms

save(Ps, file=paste(k,"Ps.Rdata",sep=""))
save(PsA, file=paste(k,"PsA.Rdata",sep=""))
save(PsAGY, file=paste(k,"PsAGY.Rdata",sep=""))
}


