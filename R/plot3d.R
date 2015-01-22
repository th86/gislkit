#Functions for three-way scatter plots

#Render a three-way scatter plot
plot3d<-function(x,y,z,colorCode=c("blue","red"),breaks=20,pch=20,cex=1,xlab="",ylab="", main=""){
    palette<-colorRampPalette(colorCode)(breaks)
    return(plot(x , y, col=palette[cut(z,breaks)], xlab=xlab,ylab=ylab,main=main, pch=pch, cex=cex))
}

#Save the plot as a TIFF file
save3d<-function(file="plot3d.tiff", x,y,z, colorCode=c("blue","red"),breaks=20,pch=20,cex=1,xlab="",ylab="",zlab="", main="", width=1024, height=1024){
    
    tiff(file=file, width=width, height=height, compression="lzw" )
    plot3d(x,y,z,colorCode,breaks,pch,cex,xlab,ylab,main)
    mtext(zlab,4)
    dev.off()

    return(0)
}
