setGeneric("drawAxis", function(GdObject, ...) standardGeneric("drawAxis"))
setMethod("drawAxis", signature(GdObject="GdObject"), function(GdObject, ...)
  return(NULL))

setMethod("drawAxis", signature(GdObject="CustomTrack"), function(GdObject,
                                                                  from,
                                                                  to, ...) {
  ylim <- displayPars(GdObject)$ylim
  hSpaceAvail <- Gviz:::vpLocation()$isize["width"]/6
  #yscale <- extendrange(r=ylim, f=0.05) #extends axis by 5%
  yscale <- ylim
  col <- Gviz:::.dpOrDefault(GdObject, "col.axis", "white")
  acex <- Gviz:::.dpOrDefault(GdObject, "cex.axis")
  acol <- Gviz:::.dpOrDefault(GdObject, "col.axis", "white")
  at <- pretty(yscale) #finds breakpoints
  #at <- at[at>=sort(ylim)[1] & at<=sort(ylim)[2]]
  if(is.null(acex))
  {
    vSpaceNeeded <- max(as.numeric(grid::convertWidth(grid::stringHeight(at),
                                                      "inches")))*length(at)*1.5
    hSpaceNeeded <- max(as.numeric(grid::convertWidth(grid::stringWidth(at),
                                                      "inches")))
    vSpaceAvail <- abs(diff(range(at)))/
      abs(diff(yscale))*Gviz:::vpLocation()$isize["height"]

    acex <- max(0.6, min(vSpaceAvail/vSpaceNeeded, hSpaceAvail/hSpaceNeeded))
  }
    nlevs <- max(1, nlevels(factor(Gviz:::.dpOrDefault(GdObject, "groups"))))
    vpTitleAxis <- grid::viewport(x=0.95, width=0.2, yscale= yscale, just=0)
    grid::pushViewport(vpTitleAxis)
    suppressWarnings(grid::grid.yaxis(gp=grid::gpar(col=acol, cex=acex), at=at))
    grid::popViewport(1)
})



##-----------------------------------------------------------------------------
## CustomTrack:
##
## A track class to allow for user-defined plotting functions
##-----------------------------------------------------------------------------
setClass("CustomTrack",
         contains=c("GdObject"),
         representation=representation(plottingFunction="function",
                                       variables="list"),
         prototype=prototype(dp=DisplayPars()))

setMethod("initialize", "CustomTrack", function(.Object, plottingFunction,
                                                variables, ...) {
  .Object <- Gviz:::.updatePars(.Object, "CustomTrack")
  .Object@plottingFunction <- plottingFunction
  .Object@variables <- variables
  .Object <- callNextMethod(.Object, ...)
  return(.Object)
})


CustomTrack <- function(plottingFunction=function(GdObject,
                                                  prepare=FALSE, ...){},
                        variables=list(), name="CustomTrack", ...){
  return(new("CustomTrack", plottingFunction=plottingFunction,
             variables=variables, name=name, ...))
}
