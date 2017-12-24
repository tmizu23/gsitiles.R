library(raster)
library(curl)
library(png)
library(jpeg)
library(spatialEco)

tileinfo<-function(x,y,zoom){
  #タイル座標からタイルの左下、右上の　WEBメルカトル座標とm/pixcelを返す
  A<-20037508.342789244
  L<- 2*A/2^zoom
  x0<-L*x-A
  y1<-A-L*y
  x1<-x0+L
  y0<-y1-L
  (list(bbox=c(x0,y0,x1,y1),res=L/256))
}

latlon2tile <- function(lon, lat, zoom) {
  #緯度、経度、ズームレベルからタイル座標を返す
  x = trunc((lon/180 + 1) * 2^zoom/2)
  y = trunc(((-log(tan((45 + lat/2) * pi/180)) + 
                pi) * 2^zoom/(2 * pi)))
  return(c(x,y,zoom))
}

latlon2tile(141.1930,42.4934,13)

loadTILE<-function(x,y,zoom,type){
  #地図画像タイルを取得
  if(type=="std"){
    baseurl<-"https://cyberjapandata.gsi.go.jp/xyz/std/"
    ext<-".png"
  }else if(type=="photo"){
    baseurl<-"https://cyberjapandata.gsi.go.jp/xyz/seamlessphoto/"
    ext<-".jpg"
  }
  myurl<-paste0(baseurl, paste(zoom, x, y, sep = "/"), ext)
  info<-tileinfo(x,y,zoom)
  print(myurl)
  if(type=="std"){
    img <- readPNG((readBin(myurl, 'raw', n=1e6)))
  }else{
    img <- readJPEG((readBin(myurl, 'raw', n=1e6)))
  }
  R<-raster(xmn=info$bbox[1], ymn=info$bbox[2],xmx=info$bbox[3], ymx=info$bbox[4], resolution=info$res,crs="+init=epsg:3857")
  G<-B<-R
  R[]<-img[,,1]*255
  G[]<-img[,,2]*255
  B[]<-img[,,3]*255
  stack(R,G,B)
}

loadDEM<-function(x,y,zoom,type){
  #png標高タイルを取得
  ext<-".png"
  if(type=="dem5a"){
    baseurl<-"https://cyberjapandata.gsi.go.jp/xyz/dem5a_png/"
  }else if(type=="dem5b"){
    baseurl<-"https://cyberjapandata.gsi.go.jp/xyz/dem5b_png/"
  }else if(type=="dem10b"){
    baseurl<-"https://cyberjapandata.gsi.go.jp/xyz/dem_png/"
  }
  myurl<-paste0(baseurl, paste(zoom, x, y, sep = "/"), ext)
  info<-tileinfo(x,y,zoom)
  print(myurl)
  img <- readPNG((readBin(myurl, 'raw', n=1e6)))
  dem<-raster(xmn=info$bbox[1], ymn=info$bbox[2],xmx=info$bbox[3], ymx=info$bbox[4], resolution=info$res,crs="+init=epsg:3857")
  dem
  dem[]<-(img[,,1]*2^16+img[,,2]*2^8+img[,,3])*255
  dem <- calc(dem, fun=function(x){ifelse(x>2^23,(x-2^24)*0.01,ifelse(x<2^23,x*0.01,x<-NA))})
}

getTile<-function(lon0,lat0,lon1,lat1,zoom,type,combine=TRUE,crop=FALSE){
  #緯度経度の左下～右上の範囲のタイルを取得
  notile=FALSE
  if(type=="dem10b" && zoom>14)
    notile=TRUE
  if((type=="dem5a"||type=="dem5b") && zoom!=15)
    notile=TRUE
  if((type=="std") && zoom>18)
    notile=TRUE
  if((type=="photo") && zoom>18)
    notile=TRUE
  if(notile){
    print("no tile for this zoom level")
    return(NULL)
  }  
  xyz0 <- latlon2tile(lon0,lat0,zoom)
  xyz1 <- latlon2tile(lon1,lat1,zoom)

  tiles<-list()
  for(y in xyz0[2]:xyz1[2]){
    for(x in xyz0[1]:xyz1[1]){
      if(type=="std"||type=="photo"){
        tile<-loadTILE(x,y,zoom,type)  
      }else{
        tile<-loadDEM(x,y,zoom,type)  
      }
      tileid<-paste(x,y,sep="-")
      tiles[[tileid]]<-tile
    }
  }

  #タイルを1つにする
  if(combine==TRUE){
    if(length(tiles)==1){
      tiles<-tiles[[1]]
    }else{ 
      names(tiles)<-NULL
      tiles$fun<-min
      tiles<-do.call(mosaic,tiles)
    }
  }
  #切り取る
  if(crop==TRUE){
    if(lon0!=lon1&&lat0!=lat1){
      latlon_r<-raster(extent(lon0,lon1,lat0,lat1),crs=CRS("+init=epsg:4612"))
      latlon_r[]<-0
      crop_r<-projectRaster(latlon_r, crs=CRS("+init=epsg:3857"))
      tiles<-crop(tiles,crop_r)
    }
  }
  return(tiles)
}

scaleAtoB<-function(x,A,B){
  #ラスタの最大値をA、最小値をBにする
  (x-minValue(x))/(maxValue(x)-minValue(x))*(A-B)+B
}

overlayColor<-function(A,B,type){
  #色を合成する
  if(type=="overlay"){
    #オーバーレイ
    f<-function(x,y){ifelse(y<128,x*y*2/255,2*(x+y-x*y/255)-255)}
  }else if(type=="multiply"){
    #乗算
    f<-function(x,y){x*y/255}
  }
  s<-c()
  for(i in 1:3){
    x<-raster(A,i)
    y<-raster(B,i)
    z<-overlay(x,y, fun=f)
    s<-c(s,z)
  }
  stack(s)
}

makeCS<-function(dem,shaded=TRUE){
  #demからCS立体図を作成して、RGBのstackを返す
  #shaded=TRUEで、陰影起伏図と合成した図を返す
  
  #曲率 白> 藍
  cur<-curvature(dem,5,type="bolstad")
  cur_r<-scaleAtoB(cur,255,0)
  cur_g<-scaleAtoB(cur,255,0)
  cur_b<-scaleAtoB(cur,255,128)
  cur_stack<-stack(cur_r,cur_g,cur_b)
  #傾斜 茶 > 白
  slp<-terrain(dem,"slope")
  slp_r<-scaleAtoB(slp,115,255)
  slp_g<-scaleAtoB(slp,66,255)
  slp_b<-scaleAtoB(slp,41,255)
  slp_stack<-stack(slp_r,slp_g,slp_b)
  #傾斜と曲率を乗算
  cs_stack<-overlayColor(slp_stack,cur_stack,type="multiply")
  if(shaded){
    #陰影起伏
    ter<-terrain(dem,c("slope","aspect"))
    hil<-hillShade(ter$slope,ter$aspect)
    hil<-scaleAtoB(hil,255,0)
    hil_stack<-stack(hil,hil,hil)
    #陰影起伏とCS立体図をオーバーレイ
    cs_stack<-overlayColor(hil_stack,cs_stack,type="overlay")
  }
  return(cs_stack)
}

