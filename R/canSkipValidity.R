setClass("canSkipValidity")
canSkipValidityInitialize <- function (.Object, ...)
{
    args <- list(...)
    skipValidity <- FALSE
    if("skipValidity" %in% names(args)) {
    	skipValidity <- args$skipValidity
	args$skipValidity <- NULL
    }
    if (length(args)) {
        Class <- class(.Object)
        ClassDef <- getClass(Class)
        snames <- allNames(args)
        which <- nzchar(snames)
        elements <- args[which]
        supers <- args[!which]
        thisExtends <- names(ClassDef@contains)
        slotDefs <- ClassDef@slots
        dataPart <- slotDefs[[".Data"]]
        if (is.null(dataPart))
            dataPart <- "missing"
        if (length(supers)) {
            for (i in rev(seq_along(supers))) {
                obj <- supers[[i]]
                Classi <- class(obj)
                if (length(Classi) > 1L)
                  Classi <- Classi[[1L]]
                if (.identC(Classi, Class))
                  .Object <- obj
                else if (extends(Classi, Class))
                  .Object <- as(obj, Class, strict = FALSE)
                else if (extends(Class, Classi))
                  as(.Object, Classi) <- obj
                else if (extends(Classi, dataPart))
                  .Object@.Data <- obj
                else {
                  extendsi <- extends(Classi)[-1L]
                  which <- match(thisExtends, extendsi)
                  which <- seq_along(which)[!is.na(which)]
                  if (length(which)) {
                    Classi <- thisExtends[which[1L]]
                    as(.Object, Classi) <- obj
                  }
                  else stop(gettextf("cannot use object of class %s in new():  class %s does not extend that class",
                    dQuote(Classi), dQuote(Class)), domain = NA)
                }
            }
        }
        if (length(elements)) {
            snames <- names(elements)
            if (anyDuplicated(snames))
                stop(gettextf("duplicated slot names: %s", paste(sQuote(snames[duplicated(snames)]),
                  collapse = ", ")), domain = NA)
            which <- match(snames, names(slotDefs))
            if (anyNA(which))
                stop(sprintf(ngettext(sum(is.na(which)), "invalid name for slot of class %s: %s",
                  "invalid names for slots of class %s: %s"),
                  dQuote(Class), paste(snames[is.na(which)],
                    collapse = ", ")), domain = NA)
            firstTime <- TRUE
            for (i in seq_along(snames)) {
                slotName <- snames[[i]]
                slotClass <- slotDefs[[slotName]]
                slotClassDef <- getClassDef(slotClass, package = ClassDef@package)
                slotVal <- elements[[i]]
                if (!.identC(class(slotVal), slotClass) && !is.null(slotClassDef)) {
                  valClass <- class(slotVal)
                  valClassDef <- getClassDef(valClass, package = ClassDef@package)
                  if (!identical(possibleExtends(valClass, slotClass,
                    valClassDef, slotClassDef), FALSE))
                    slotVal <- as(slotVal, slotClass, strict = FALSE)
                }
                if (firstTime) {
                  slot(.Object, slotName, check = FALSE) <- slotVal
                  firstTime <- FALSE
                }
                else {
                  `slot<-`(.Object, slotName, check = FALSE,
                    slotVal)
                }
            }
        }
        if(!skipValidity) validObject(.Object)
    }
    .Object
}
#' @export
#' @rdname initialize
setMethod(f = "initialize", signature = "canSkipValidity", definition = canSkipValidityInitialize)
