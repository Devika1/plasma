readIn<-function(nme){
  if(!file.exists(nme)) stop(paste("file not exists", nme))
  v = fromJSON(nme)
  x = as.numeric(names(v)) 
  t(rbind(x,v))
}

makeCumul<-function(res){
  di = dim(res)
  for(i in 2:(di[1])){
    res[i,2] = res[i-1,2] + res[i,2] 
  }
  res
}

do_smooth<-function(clusters, xlim = c(2,2000), inds = 2:(dim(clusters)[2]),span = 0.05){
  if(span==0) return (clusters)
  nmes = names(clusters)
  #print(nmes)
  #print(dim(clusters))
  for(j in inds){	
    names(clusters)[j] = "depth"
    
    if(max(clusters$depth,na.rm=T)>0){
      # y = try(.loess(clusters, span))
      y = ma(clusters[,j],span)
      # try(.smooth(clusters[,c(1,j)], xlim=NULL,span,F))
      if(inherits(y,"try-error")) {
        print("error loess")
        #print(which(is.na(clusters$pos)))
        #plot(clusters$pos,clusters$depth)
        #print()
      }else{
        #print("is ok ")
        y[clusters$depth==0] =0
        clusters[,j]=  y
      }
    } 
    names(clusters)[j] = nmes[j]
  }
  
  clusters
}

.loess<-function(clusters, span){
  cars.lo <- stats::loess(depth ~ x, clusters, span = span)
  y = try(stats::predict(cars.lo, clusters, se = FALSE))
  y
}

calcCell<-function(a1){
  num = (exp(a1$filtered$merged) - exp(a1$filtered$normal))
  denom = (exp(a1$filtered$tumor) - exp(a1$filtered$normal))
  res = num/denom
  res[which(num<0)] = NA
  res[which(denom<0)] = NA
  res
}


countLt<-function(res, thresh, lt = T, 
                  secondThresh1=if(lt) rep(min(thresh)-1, length(thresh)) else rep(max(thresh)+1, length(thresh)), 
                  exclude = c()
){
  secondThresh = secondThresh1;
  if(length(secondThresh1)==1) secondThresh = rep(secondThresh, length(thresh))
  res1 = rep(0,length(thresh))
  res[which(thresh %in% exclude),1] = 0
  for(i in 1:length(thresh)){
    if(lt) res1[i] = sum(res[which(res[,1] <thresh[i] & res[,1]>= secondThresh[i]),2])
    else res1[i] = sum(res[which(res[,1] >=thresh[i] & res[,1]<secondThresh[i]),2])
  }
  res1
}
#apply(l$c1524s$c[c(100,200),],2,normEnd)[1,]

#w = c((proportion of non tumour reads < 100bp) ~0.003  , proportion of tumour reads  < 100bp) ~0.079
#w = c((proportion of non tumour reads < 100bp) ~0.003  , proportion of tumour reads  < 50bp) ~0.003275935
calcC1<-function(c, w){
  num = w[2] * c
  denom = (w[2] - w[1])*c + w[1]
  num/denom 
}

calcC<-function(c1,w){
  num = w[1] *c1
  denom = w[2] + (w[1] - w[2])*c1
  num/denom
}

#b is prob of mutated read in less than 100bp reads
#a is prob of mutated read overall
#w = c((proportion of non tumour reads < 100bp) ~0.003  , proportion of tumour reads  < 100bp) ~0.079
#rates = c(prob mutated read < 100bp, prob mutated read overall)
#mutrates = c(prob mutated read in tumour sample, prob mutated read in normal sample)
#K indicates increase in prob of selecting mutated read in longer reads (ie in <200bp vs < 100bp)
getC<-function(w,rates, mutrates = c( 0.1, 1e-5), K = 2){
  x = (1:2000)/10000  # up to 20%
  y = calcC1(x,w)
  y1  = (y-x)
  target = (rates[1]-rates[2]/K)/(mutrates[1]-mutrates[2])
  y2 = abs(y1 - target)
  x[which(y2==min(y2))]
}





normEnd<-function(v){
  v/v[length(v)] 
}
ratio<-function(v){
  v[2]/v[1]
}

analy<-function(l, col, c=0.01, points= 80:140){
  comb = l[[col]]$c
  y = rep(0, length(points))
  y1 = rep(0,length(points))
  for(i in 1:length(points)){
    w100 = apply(comb[c(points[i],200),],2,normEnd)[1,]
    c1 = calcC1(0.01, w100)
    y1[i] = c1*100;
    y[i] = w100[1]*100;
  }
  
  if(col==1) plot(c(points, points), c(y,y1),col='white');
  lines(points, y, col = col, lty = 1)
  lines(points,y1, col  = col, lty = 2)
  if(col==3){
    legend(130, y =15, c('1524', '1249', '1494'), text.col = 1:3)
  }
}
.avgInner<-function(v){
  c(min(v), mean(v), max(v))
}

readDataAll<-function(names, ext, x = 2:1999,  nme1="q",probs = c(0.0,0.5, 1.0)){
  res = array(0,dim = c(length(x), length(names)), dimnames = list(x,names))
  for(i in 1:length(names)){
    inn = names[i]
    nme_all =paste0(inn,ext)
    r_i  = readIn(nme_all);
    ind1 = r_i[,1] %in% x
    
    res[ match(r_i[ind1,1], x),i] = r_i[ind1,2]
  }
  combined = apply(res,1,sum)
  df= data.frame(res)
  df
}
.addAvg<-function(ratios, txt,probs = c(0,0.5,1.0)){
  #txt = "less than"
  
  r1 = ratios[,grep(txt, names(ratios))]
  df1 = t(apply(r1,1,quantile,probs= probs,na.rm=T))
  dimnames(df1)[[2]] = paste(txt,"q",probs*100,sep="_")
  data.frame(cbind(r1,df1))
}
ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
makeDensity<-function(vec, smooth){
  r = vec/(sum(vec))
  if(smooth>0) r = ma(r,smooth)
  r
}

ident<-function(x) x
cumsum1<-function(x) 1- cumsum(x)
plotAll1<-function(x,li, n = 5, log = "y", inds = 1:length(x), func =ident){
  
  colvec = c()
  linevec = c()
  nmevec = c()
  y = func(ma(li[[1]][inds,1],n=n))
  inds1 = which(!is.na(y))
  maxy = max(y, na.rm=T)
  firsty = y[inds1[1]]
  main  = paste("Density plot (",n,"bp sliding window)", sep="" ) 
  if(maxy>0.99){
    
    if(firsty<0.001) main  = "Cumulative density plot" else main = "1 - Cumulative density plot"
  }
  
  plot(x[inds][inds1], y[inds1], type='l', col= 0, lty=1, log=log, xlab="Fragment length(bp)", ylab="Density", main = main)
  for(i in 1:length(li)){
    li1 = li[[i]]
    for(j in 1:dim(li1)[[2]]){
      lines(x[inds], func(ma(li1[inds,j],n=n)), type='l', col= i, lty=j)
      colvec = c(colvec,j)
      linevec = c(linevec, i)
      nmevec = c(nmevec, paste(names(li)[i], dimnames(li1)[[2]][j]))
    }
  }
  maxx = max(x[inds][inds1])*0.7
  if(log=="xy") maxx = min(x[inds][inds1])
  print(paste(maxx, maxy))
  
  legend(maxx, y = maxy, nmevec,lty = colvec, col = linevec, cex=0.7)
  
}
readData<-function(inn = "", shared = FALSE, show = T, lines=F, col=1, 
                   ext_all  =  '_P0_read_length_count.json',lt=T, yback = NULL,
                   #ext_unique = if(shared) '_P0_shared_read_length_count.json' else '_P0_merged_read_length_count.json',	
                   ext_unique = if(shared) '_P0_shared_read_length_count.json' else '_P0_merged_read_length_count.json', 
                   ylim = NULL,		
                   
                   #ylim = if(shared) c(0,0.2) else c(0,3),
                   targetinds = c(100,200),
                   thresh1 = 1:2000, 
                   secondThresh=if(lt) 0 else 20000,
                   norm = NULL , lty = 1,  main = ""
){
  
  if(lt) {
    'Proportion of reads with mutation in fragments < x'
    if( secondThresh[1] >0) main = paste(main, 'and fragment >', secondThresh);
  }else{
    main ='Proportion of reads with mutations in fragments >x'
    if(secondThresh[1] <2000) main = paste(main, 'and fragment <=', secondThresh);
  }
  if(length(secondThresh)>1){
    window = thresh[1] - secondThresh[1]-1
    main = paste0('Proportion of reads with mutation in fragment windows of size ', window)
    
  }else{
    window = 0;
  }
  nme_all =paste0(inn,ext_all)
  nme_unique = paste0(inn,ext_unique) 
  print(nme_all)
  print(nme_unique)
  all  = readIn(nme_all);
  unique= readIn(nme_unique);
  allcount = countLt(all, thresh1, lt = lt, secondThresh = secondThresh)
  uniquecount = countLt(unique, thresh1, lt=lt, secondThresh = secondThresh)
  combined = cbind(allcount, uniquecount)
  dimnames(combined) = list(thresh1, c(nme_all, nme_unique))
  diff =  apply(log(combined),1,diff)
  y = exp(diff)
  ylab = 'proportion  of reads with mutation';
  if(!is.null(norm)) ylab = 'normalised prlLT = plotAll(outf=NULL, lt=T, secondThresh = 0, thresh=0:300, norm=NULL, ylim = c(0,2.2))oportion (%) of reads with mutation';
  if(!is.null(norm)) {
    y = y/y[if(lt) length(thresh1) else 1]
  }
  if(is.null(ylim)) ylim = c(min(y[y>0], na.rm=T)/10, max(y, na.rm=T) * 1.5)
  if(ylim[1]<1e-20) ylim[1] = 1e-20
  if(show & length(thresh1)>2) {
    if(lines){
      lines (thresh1, y, xlab = 'threshold', ylab = ylab, main = main, type = 'l', col=col, lty = lty)
      if(!is.null(yback)){
        ratios = apply(cbind(yback,y),1,ratio)
        #print(ratio)
        lines (thresh1, ratios, xlab = 'threshold', ylab = ylab, main = main, type = 'l', col=col, lty = lty+1)
      }
    }else{
      plot(thresh1, y, xlab = 'threshold', ylab = ylab,   main = main, type = 'l', col=col, ylim = ylim, lty = lty, log='y')
    }
  }
  inds = match(targetinds, thresh1)
  w =  apply(combined[inds,],2,normEnd)[1,]
  rates = apply(combined[inds,],1,ratio)
  list(combined=combined, w = w, thresh = thresh1, rates = rates,all = all, allcount = allcount, uniquecount = uniquecount, unique = unique, nme_all = nme_all, nme_unique = nme_unique, y = y)
}
calcDiff<-function(combined){
  diff =  apply(log(combined),1,diff)
  print(diff[which(thresh==200)] - diff[which(thresh==100)] )
  cbind(thresh1, diff)
}

average<-function(l, 
                  groups = list("merged"= c("c1524", "c1249", "c1494"),"tumor"= c("c1524s", "c1249s", "c1494s"),  "normal"= c("c609", "c616"))
){
  thresh = l[[1]]$thresh
  res = list()
  leng = length(l[[1]]$uniquecount)
  for( i in 1:length(groups)){
    inds = which(names(l) %in% groups[[i]])
    uniq = rep(0, leng)
    all =  rep(0, leng)
    for(j in 1:length(inds)){
      uniq = uniq + l[[inds[j]]]$uniquecount
      all = all + l[[inds[j]]]$allcount
    }
    res[[i]] = list(uniquecount = uniq, allcount = all, thresh = thresh)
  }
  names(res) = names(groups)
  res
}
plotlines<-function(l, xlim = NULL, outf=NULL, ylim = NULL){
  if(!is.null(outf))pdf(outf)
  if(is.null(xlim)) xlim = c(min(l[[1]]$thresh),max(l[[1]]$thresh))
  filtered= list()
  all = list()
  thresh = l[[1]]$thresh
  y_f   = log(l[[1]]$uniquecount/max(l[[1]]$uniquecount));
  y_all = log(l[[1]]$allcount/max(l[[1]]$allcount))
  plot(thresh, y_all, type='l', col=1, lty=1, xlim = xlim, ylim = ylim)
  lines(thresh, y_f, type='l', col=1, lty=2)
  all[[1]] = y_all
  filtered[[1]] = y_f
  for(i in 2:length(names(l))){
    y_all = log(l[[i]]$allcount/max(l[[i]]$allcount))
    y_f = log(l[[i]]$uniquecount/max(l[[i]]$uniquecount))
    all[[i]] = y_all
    filtered[[i]] = y_f
    lines(thresh, y_all, type='l', col=i, lty=1)
    lines(thresh, y_f, type='l', col=i, lty=2)
    
  }
  names(all) = names(l)
  names(filtered) = names(l)
  maxy = log(0.9);
  if(!is.null(ylim)) may=ylim[2]
  legend(xlim[2]*0.7, y = maxy, names(l), text.col = 1:(length(names(l))))
  if(!is.null(outf)) dev.off()
  list(all=all, filtered = filtered)
}

plotDiff<-function(a1, thresh, outf = NULL){
  if(!is.null(outf))pdf(outf)
  diff = list()
  diff[[1]] = a1$filtered$tumor - a1$filtered$normal
  diff[[2]] = a1$all$tumor - a1$all$normal
  diff[[3]] = a1$filtered$tumor - a1$filtered$merged
  names(diff) =c( "filtered tumor vs filtered normal","all_tumor vs all_normal",  "filtered tumor vs filtered merged")
  plot(thresh, diff[[1]], type = 'l', col=1, lty = 1)
  lines(thresh, diff[[2]], type = 'l', col = 2, lty = 1)
  lines(thresh, diff[[3]], type = 'l', col = 3, lty = 1)
  legend(thresh[length(thresh)*0.5] , y = max(diff[[1]], na.rm=T), names(diff), text.col = 1:3)
  if(!is.null(outf))dev.off()
}

plotAll<-function(norm=NULL, thresh = 1:2000, outf = NULL, lt=T, ylim = NULL,
                  targetinds = c(100,200),
                  secondThresh=if(lt) 0 else 20000
){
  if(!is.null(outf))pdf(outf, width = 21, height = 7)
  shared = F
  l = list();
  names = c("1524", "1249", "1494", "609", "616");
  l[[1]] = readData(names[1], thresh = thresh, shared = shared, norm = norm, lt = lt, lines=F, col = 1, ylim = ylim, secondThresh=secondThresh, targetinds = targetinds)
  for(i in 2:length(names)){
    l[[i]] = readData(names[i], thresh = thresh, shared = shared, norm = norm, lt = lt, lines=T, col = i, secondThresh=secondThresh, targetinds = targetinds)
  }
  shared = T
  leng = length(names);
  for(i in 1:3){
    l[[i+ leng]]  = readData(names[i], thresh = thresh, lines=T,lty=2, col=i, shared = shared, norm = norm, lt = lt, secondThresh=secondThresh, targetinds = targetinds, yback = l[[i]]$y)
  }
  legend(150, y = max(l[[1]]$y, na.rm=T), c('1524', '1249', '1494', '065'), text.col = 1:5)
  if(!is.null(outf)) dev.off()
  names(l)  = c(paste0('c', names), paste0('c', names[1:3], 's'))
  l
}
.smooth<-function(d1, xlim = NULL, step = 10, cumul = F){
  
  if(!is.null(xlim)){
    inds1 =  which(d1[,1]>=xlim[1] & d1[,1] <=xlim[2])
    d = d1[inds1,,drop=F]
    todo = xlim[1]:xlim[2]
  }else{
    todo = 2:1998
  }
  d = cbind(todo, rep(0, length(todo)))
  mi = match(d1[,1], todo)
  d[mi[!is.na(mi)],2] = d1[!is.na(mi),2]
  if(step>10){
    max = dim(d)[1]
    vec = seq(1, max, step)
    matr = array(NA,dim = c(length(vec), 2))
    for(i in 1:length(vec)){
      inds = vec[i]:min(vec[i]+step, max)
      matr[i,] =  apply(d[inds,, drop=F],2,sum)
      matr[i,1] = matr[i,1]/length(inds)
    }
  }else{
    matr = d
  }
  if(cumul) matr[,2] = .cumul(matr[,2])
  matr[,2]
}
ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}

.cumul<-function(vec) {
  for(i in 2:length(vec)) vec[i] = vec[i] + vec[i-1];
  vec
}


replot<-function(l,  xlim = c(100,200), nme = names(l)[1:5], norm=T, samesum = F, smooth = 1, cumul = F, log="y", normtop = F, ymin = 1, secondnorm=F){
  
  plott = F
  #plot(x,y, log=log, type='l',ylim = ylim,  xlim = xlim, col="white")
  for(i in 1:length(nme)){
    inds = grep(nme[i], names(l))
    lt = 1
    print(inds)
    ytop = NULL
    for(j in 1:length(inds)){
      if(j==1){
        d = .smooth(l[[inds[j]]]$all, xlim, smooth ,cumul)
        x = d[,1]
        y =  d[,2]
        if(j==1) ytop = y
        if(!samesum || j==1) sumy = sum(y)
        if(normtop) y = y/ytop
        else if(norm) y = y/sumy
        if(plott==F){
          miny = min(ymin,min(y))
          if(norm) miny =min(ymin/sumy,min(y)) 
          ylim = c(miny,max(y))
          plot(x,y, log=log, type='l',ylim = ylim,  xlim = xlim, col=i, lty=lt)
          plott=T
        }
        lines(x,y,   type='l', col=i, lty=lt)
        lt = lt+1
        ytop =  d[,2]
      }
      d = .smooth(l[[inds[j]]]$unique, xlim, smooth, cumul)
      x = d[,1]
      y =  d[,2]
      print(length(ytop)) 
      print(length(y))
      if(!samesum) sumy = sum(y)
      if(normtop) y =y/ytop
      else	if(norm) y = y/sumy
      
      lines(x,y,  type='l', lty=lt, col=i)
      lt = lt +1
      if(secondnorm) ytop = d[,2]
    }
  }
  legend(xlim[1], y =ylim[2], nme, text.col = 1:length(nme))
  #lines(l[[1]]$unique,  type='l', col=3)
}


