source("gsitiles.R")

#-------------------------
#おっぱい山の可視化
#-------------------------

#おっぱい山データベースの読み込み
oppai_csv<-read.csv("https://gist.githubusercontent.com/tmizu23/29c24df1662bbf976e2407f843f43c15/raw/1c0d7cf03e2e44ecad41d7b21772820e66bc95b3/oppaiyama.csv",fileEncoding="UTF-8")
offsets<-c(0.009,0.005,0.005,0.005,0.004,0.004,0.005,0.005,0.005,0.005)
stdzoom<-c(14,15,15,15,17,17,15,15,15,15)

for(i in 1:nrow(oppai_csv)){
  #おっぱい山の数だけCS立体図と地図画像をPNGで書き出し
  offset<-offsets[i]
  opp<-oppai_csv[i,]
  zoom<-stdzoom[i]
  dem<-getTile(opp$lon-offset,opp$lat-offset,opp$lon+offset,opp$lat+offset,14,type="dem10b",crop=TRUE)
  cs<-makeCS(dem,shaded=TRUE)
  stdmap<-getTile(opp$lon-offset,opp$lat-offset,opp$lon+offset,opp$lat+offset,zoom,type="std",crop=TRUE)
  png(paste("oppaiyama",i,".png",sep=""),600,340)
  par(mfrow=c(1,2),oma=c(0,0,2,0))
  plotRGB(cs)
  plotRGB(stdmap)
  titlestr<-paste(oppai_csv$name[i],"(",oppai_csv$area[i],")",sep="")
  title(titlestr,outer=T, cex=2)
  dev.off()
}

#----------------------
#その他、使用例
#----------------------

#標高データを取得・表示・投影変換・書き出し
dem<-getTile(137.665215,36.294582,137.7,36.314582,15,type="dem5a")
plot(dem)
demUTM53<-projectRaster(dem, crs=CRS("+init=epsg:3099"))
writeRaster(demUTM53,"demUTM53.tif")

#陰影起伏図の作成・表示
slope<-terrain(dem,"slope")
aspect<-terrain(dem,"aspect")
shade<-hillShade(slope,aspect)

plot(shade,col=grey(0:100/100),legend=F)

#標準地図とCS立体図を乗算合成・表示・geoTIFF書き出し
cs<-makeCS(dem,shaded=TRUE)
stdmap<-getTile(137.665215,36.294582,137.7,36.314582,15,type="std")
stdmap_cs<-overlayColor(stdmap,cs,type="multiply")
plotRGB(stdmap_cs)
writeRaster(stdmap_cs,"stdmap_cs.tif",datatype="INT1U") #datatypeを指定

#写真と陰影起伏を透過表示
photo<-getTile(137.665215,36.294582,137.7,36.314582,15,type="photo")
plot(shade,col=grey(0:100/100),legend=F)
plotRGB(photo,alpha=180,add=T)

#写真と陰影起伏を乗算表示
shade2<-scaleAtoB(shade,255,0)#RGB表示用にスケール変換
shade_stack<-stack(shade2,shade2,shade2)#RGB用にstack
photo_shade<-overlayColor(shade_stack,photo,type="multiply")
plotRGB(photo_shade)
writeRaster(photo_shade,"photo_shade.tif",datatype="INT1U") #datatypeを指定
