\name{SpectralMap}
\alias{SpectralMap}
\title{
Diffusion Map and Spectral Map
}
\description{
Implements the diffusion map method of dimensionality reduction and spectral method of combining multiple diffusion maps, including creation of the spectra and visualization of maps.
}
\usage{
SpectralMap(data, epsilon=0.1, range=1, Plot2D=FALSE, Plot3D=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
  	a list of datasets and each column in each dataset is a variable
}
  \item{epsilon}{
  	parameter in the Gaussian kernel
}
  \item{range}{
  	indexes of the datasets in the data list to be combined and computed. If length(range)==1, only diffusion map will be computed. Otherwise, spectral map will be computed
}
  \item{Plot2D}{
		a logical value indicating whether a 2-D map should be plotted
}
  \item{Plot3D}{
		a logical value indicating whether a 3-D map should be plotted
}
}
\value{
\item{singularvector}{
		the spectra of either diffusion map or spectral map
}
}
\examples{
#generate two datasets
n <- 100
theta <- 2*pi*seq(from=0, to=1-1/n, by=1/n)
r = (1 + cos(6*theta)/4)

# X is a circle
X1 = cos(theta)
X2 = sin(theta)
X<-data.frame(X1,X2)

#Y is a hexagon
Y1 = r*cos(theta)
Y2 = r*sin(theta)
Y<-data.frame(Y1,Y2)

#create data list
Data<-list(X,Y)

#create the diffusion map of X
example1<-SpectralMap(Data, epsilon=0.1, range=1, Plot2D=TRUE, Plot3D=FALSE)
#create the spectral map from X to Y
example2<-SpectralMap(Data, epsilon=0.1, range=1:2, Plot2D=TRUE, Plot3D=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.