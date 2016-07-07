SpectralMap<-function(data, epsilon=0.1, range=1, Plot2D=FALSE, Plot3D=FALSE)
{
  if(length(range)==1){
    # Pairwise Euclidean Distance
    # D takes values 0 to infty
    # Wehre 0 mean very similar infty means very different
    D<-rdist(data[[range]],data[[range]])
    # Affinity
    # K takes values 0 to 1
    # where 1 means very similar 0 means very different
    # K is the gaussian kernal
    K <- exp(-D^2/epsilon)
    d <- rowSums(K)
    D <- diag(1/d)
    P <- D %*% K
    
    # SVD decomposition
    U <- svd(P)$u
  }
  if(length(range)>=1){
    # Pairwise Euclidean Distance List
    # D takes values 0 to infty
    # Wehre 0 mean very similar infty means very different
    D<-list()
    for(r in range[1]:range[length(range)])
    {
    D[[r]]<-rdist(data[[r]],data[[r]])
    }
    # Affinity
    # K takes values 0 to 1
    # where 1 means very similar 0 means very different
    # K is the gaussian kernal
    K<-list()
    d<-list()
    for(i in range[1]:range[length(range)])
    {
      K[[i]] <- exp(-D[[i]]^2/epsilon)
      d[[i]] <- rowSums(K[[i]])
    }
    for(i in range[1]:range[length(range)])
    {
      D[[i]] <- diag(1/d[[range[1]]])
    }
    P <- list()
    for(i in range[1]:range[length(range)])
    {
      P[[i]]<- D[[i]] %*% K[[i]]
    }
    # Spectral matrix is defined as matrix chain multiplication
    SpectralP<-diag(1,dim(data[[r]])[1])
    for(i in range[1]:range[length(range)])
    {
      SpectralP<- SpectralP %*% P[[i]]
    }
    # SVD decomposition
    U <- svd(SpectralP)$u
  }
  
# 2D-plot
if(Plot2D){plot(U[,2],U[,3],xlab="",ylab="")}
# 3D-plot
if(Plot3D){scatterplot3d(U[,2],U[,3],U[,4],xlab="",ylab="",zlab="")}
# Return Singular Value
object<-list(singularvector=U)
object
}