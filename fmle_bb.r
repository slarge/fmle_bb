# fmle_bb()    {{{

setGeneric('fmle_bb', function(object, start, ...)
  standardGeneric('fmle_bb'))


setMethod('fmle_bb',
          signature(object='FLModel', start='FLPar'),
          function(object, start, method='Nelder-Mead', fixed=list(),
                   control=list(trace=1), lower=rep(-Inf, dim(params(object))[2]),
                   upper=rep(Inf, dim(params(object))[2]), ...)
          {
            values <- as.list(FLCore::iter(start,1))
            names(values) <- dimnames(start)$params
            
            fmle_bb(object, values, method, fixed, control, lower, upper, ...)
          }
)

# fmle {{{
setMethod("fmle_bb", signature(object="FLSR", start="ANY"),
          function(object, start, ...)
          {
            res <- callNextMethod()
            # AR1 models
            if('rho' %in% dimnames(params(object))$params)
            {
              n <- dim(rec(res))[2]
              rho <- c(params(res)['rho',])
              residuals(res) <- as.numeric(NA)
              residuals(res)[,-1] <- (rec(res)[,-1] - rho*rec(res)[,-n] - fitted(res)[,-1] +
                                        rho*fitted(res)[,-n])
            }
            # lognormal models
            else if(object@logerror)
              residuals(res) <- log(rec(res)) - log(fitted(res))
            return(res)
          }
) # }}}


setMethod('fmle_bb', 
          signature(object='FLModel', start='ANY'),
          function(object, start, method='Nelder-Mead', fixed=list(),
                   control=list(trace=1), lower=rep(-Inf, dim(params(object))[1]),
                   upper=rep(Inf, dim(params(object))[1]), seq.iter=TRUE, preconvert = FALSE, 
                   inParallel = TRUE, ...)
          { 
            
            ## Figure out what should be pre-processed locally and how to export to BB ##
            # rm(list = ls())
            # data(ple4)
            # ple4SR<-as.FLSR(ple4)
            # #### Specifying the stock recruitment relationship and error model
            # model(ple4SR)<-bevholt()
            # ple4SR <- propagate(ple4SR, iter = 1000)
            # object <- ple4SR
            # method='Nelder-Mead'
            # fixed=list()
            # seq.iter=TRUE
            # control=list(trace=1)
            # library(doParallel)
            # library(FLCore)
            # #
            
            # TODO Check with FL
            args <- list(...)
            call <- sys.call(1)
            logl <- object@logl
            
            # get parameter names by matching elements in param slot
            parnm <- names(formals(logl))[names(formals(logl))%in%
                                            dimnames(object@params)$param]
            
            # get fixed parameter names
            fixnm <- names(fixed)
            # fixed must match params
            if(any(!fixnm %in% parnm)) {
              stop("some named arguments in 'fixed' are not arguments to the
                   supplied log-likelihood")
            }
            # HACK! clean up fixed list if elements are named vectors
            fixed <- lapply(fixed, function(x){ names(x) <- NULL; x})
            
            # create list of input data
            #   get FLQuant slots' names
            datanm <- getSlotNamesClass(object, 'FLArray')
            # Include FLQuants contents too
            flqs <- getSlotNamesClass(object, 'FLQuants')
            for (i in length(flqs)) {
              datanm <- c(datanm, names(slot(object, flqs[i])))
            }
            datanm <- c(datanm, getSlotNamesClass(object, 'numeric'))
            #   get those in formals of logl
            datanm <- datanm[datanm%in%names(formals(logl))]
            
            # limits
            if(method %in% c('L-BFGS-B', 'Brent'))
            {
              if(missing(lower) && !is.null(lower(object)))
                # if is(lower, function)
                lower <- lower(object)[match(parnm, names(fixed), nomatch=0)==0]
              if(missing(upper) && !is.null(upper(object)))
                upper <- upper(object)[match(parnm, names(fixed), nomatch=0)==0]
            } else
            {
              lower <- -Inf
              upper <- Inf
            }
            
            # gr function
            if(!is.null(body(object@gr)))
            {
              gr <- function(par)
              {
                pars <- as.list(par)
                names(pars) <- names(start)
                pars[fixnm] <- lapply(fixed, FLCore::iter, it)
                return(-1*(do.call(object@gr, args=c(pars, data))))
              }
            } else
              gr <- NULL
            
            # create logl function
            loglfoo <- function(par) {
              pars <- as.list(par)
              names(pars) <- names(start)
              pars[fixnm] <- lapply(fixed, FLCore::iter, it)
              return(-1*(do.call(logl, args=c(pars, data))))
            }
            
            # input data
            alldata <- list()
            # slots
            for(i in datanm[!datanm %in% names(covar(object))]) {
              alldata[[i]] <- slot(object, i)
            }
            
            if(length(covar(object)) > 0) {
              for (i in datanm[datanm%in%names(covar(object))]) {
                alldata[[i]]  <- covar(object)[[i]] 
              }
            }
            
            # add dimnames if used
            dimna <- dimnames(slot(object, datanm[1]))[names(slot(object, datanm[1]))%in%
                                                         all.vars(object@model)]
            if(length(dimna) > 0)
            {
              # get them in the right shape
              dimdat <- lapply(dimna, function(x)
              {
                out <- slot(object, datanm[1])
                out[] <- as.numeric(x)
                return(out)
              })
              alldata <- c(alldata, dimdat)
            }
            
            # iterations
            if(seq.iter) 
            {
              iterReps <- dims(object)$iter
              # iters in fixed
              if(length(fixnm) >= 1)
              {
                fiter <- unlist(lapply(fixed, length))
                if(!all(fiter == 1))
                {
                  fiter <- fiter[fiter > 1]
                  # all iters in fixed are equal?
                  if(any(fiter/fiter[1] != 1))
                    stop("objects in fixed have different number of iters")
                  # are iter in object 1 and fixiter > 1? use fixiter
                  if(iterReps == 1 & fiter > 1)
                    iterReps <- fiter
                  # are they different and > 1? STOP
                  else if(fiter > 1 & fiter != iterReps)
                    stop("different iters in fixed and object")
                }
              }
            } else { 
              iterReps <- 1
            }
            
            # logLik
            logLik <- rep(NA, iterReps)
            class(logLik) <- 'logLik'
            attr(logLik, 'df') <- length(parnm) - length(fixed)
            object@logLik <- logLik
            
            # Correct FLPar, fitted and residuals
            if(iterReps > dim(object@params)[length(dim(object@params))])
            {
              params(object) <- FLPar(iter=iterReps, params=dimnames(object@params)$params)
            }
            
            fitted(object) <- propagate(fitted(object), iterReps)
            residuals(object) <- propagate(residuals(object), iterReps)
            
            # vcov
            object@vcov <- array(NA, dim=c(rep(length(parnm)-length(fixed),2), iterReps),
                                 dimnames=list(parnm[!parnm%in%names(fixed)],parnm[!parnm%in%names(fixed)],
                                               iter=1:iterReps))
            object@hessian <- object@vcov
            
            # for (it in 1:iterReps) {
            # data
            # if(seq.iter){
            #   data <- lapply(alldata, FLCore::iter, it)
            # } else {
            #   data <- alldata
            # }
            
            # do preconversion of data objects
            # if(preconvert) {
            #   data <- lapply(data, c)
            # }
            # 
            # start values
            if(missing(start)) {
              # add call to @initial
              if(is.function(object@initial)) {
                start <- lapply(1:iterReps, function(x) as(do.call(object@initial, 
                                                                   args = lapply(alldata, 
                                                                                 FLCore::iter, 
                                                                                 x)[names(formals(object@initial))]),
                                                           'list'))

                # start <- as(do.call(object@initial, 
                #                     args = data[names(formals(object@initial))]),
                #             'list')
              } else {
                start <- formals(logl)[names(formals(logl))%in%parnm] 
              }
            } else {
              # HACK! clean up fixed list if elements are named vectors
              start <- lapply(start, function(x){ names(x) <- NULL; x})
            }
            # MAKE SURE ThiS STILL WORKS
            if(!is.null(fixnm)){
              start[fixnm] <- NULL
            }
            
            if(any(!names(start) %in% parnm)) {
              stop("some named arguments in 'start' are not arguments to the
                     supplied log-likelihood")
            }
            
            ## START PARALLEL STUFF ##
            unregister <- function() {
              env <- foreach:::.foreachGlobals
              rm(list=ls(name=env), pos=env)
            } # close unregister function

            if(inParallel == TRUE) {
              #
              detectedCores <- parallel::detectCores() - 2
              cl <- parallel::makeCluster(detectedCores)
              doParallel::registerDoParallel(cores = cl)
              #
            } # Close inParallel == TRUE
            # 
            # Register a sequential backend
            if(inParallel == FALSE) {
              foreach::registerDoSEQ()
            } # Close inParallel == FALSE
            # 
            
            # it = 1
            # dat <- list(start = start,
            #             alldata = alldata,
            #             loglfoo = loglfoo,
            #             parnm = parnm,
            #             gr = gr)
            
            out <- foreach(it = 1:iterReps,
                           # .export = c("alldata", "start", "parnm", "loglfoo", "gr"),
                           .packages = c("FLCore")) %dopar% {

                             data <- lapply(alldata, FLCore::iter, it)
                             startit <- start[[it]][order(match(names(start[[it]]), parnm))]

                             # add small number to start if 0
                             startit <- lapply(startit, function(x) if(x == 0) x/100000 else x)

                             if(is.null(startit)) {
                               stop("No starting values provided and no initial function available")
                             }

                             # TODO protect environment
                             out <- do.call('optim', c(list(par = unlist(startit),
                                                            fn=loglfoo,
                                                            # parnm = parnm,
                                                            method=method,
                                                            hessian=TRUE,
                                                            control=control,
                                                            lower=lower,
                                                            upper=upper,
                                                            gr=gr)))
                             out$iter <- it

                             # warning if convergence is not 0, and do not load results
                             if(out$convergence != 0) {
                               warning("optimizer could not achieve convergence")
                             }
                             return(out)
                           } # Close FOREACH
            
            if(inParallel == TRUE){
              parallel::stopCluster(cl = cl)
            } # close inParallel == TRUE
            if(inParallel == FALSE){ 
              foreach::registerDoSEQ()
            } # close inParallel == FALSE
            unregister()
            stopImplicitCluster()
            # 
            
            # output
            for(ir in 1:iterReps) {

              # place out$par in right iter dim
              FLCore::iter(object@params[names(start[[ir]]),], ir) <- out[[ir]]$par

              # fixed
              if(length(fixed) > 0) {
                FLCore::iter(object@params, ir)[fixnm,] <- unlist(lapply(fixed,
                                                                         FLCore::iter,
                                                                         ir))
              }

              # TODO make details list of lists if iter > 1?
              # FLCore::iter(object@details, ir) <- list(call="call",
              #                                          value=out[[ir]]$value,
              #                                          count=out[[ir]]$counts,
              #                                          convergence=out[[ir]]$convergence,
              #                                          message=out[[ir]]$message)

              object@details <- list(call=call, value=out[[ir]]$value,
                                     count=out[[ir]]$counts,
                                     convergence=out[[ir]]$convergence,
                                     message=out[[ir]]$message)

              # vcov & hessian
              coef <- out[[ir]]$par
              object@vcov[,,ir] <-
                if (length(coef))
                {
                  if(det(out[[ir]]$hessian) != 0)
                  {
                    tmphess <- try(solve(out[[ir]]$hessian), silent=TRUE)
                    if(class(tmphess) =='try-error')
                    {
                      matrix(numeric(0), length(coef), length(coef), dimnames=list(names(coef),
                                                                                   names(coef)))
                    } else
                      tmphess
                  } else
                    0
                } else
                  0
              object@hessian[,,ir] <- -out[[ir]]$hessian

              # logLik
              object@logLik[ir] <- -out[[ir]]$value
              attr(object@logLik, 'nobs') <- length(lapply(alldata, FLCore::iter, ir)[[1]])

              # fitted & residuals
              FLCore::iter(fitted(object), ir) <- predict(FLCore::iter(object, ir))
              FLCore::iter(residuals(object), ir) <- FLCore::iter(slot(object,
                                                                       as.list(object@model)[[2]]), ir) - FLCore::iter(fitted(object), ir)


              # force dimnames[1:5] in 'fitted' and 'residuals' to match
              dimnames(fitted(object))[1:5] <- dimnames(do.call(as.character(as.list(object@model)[2]),
                                                                list(object)))[1:5]
              dimnames(residuals(object)) <- dimnames(fitted(object))

            } # CLOSE output loop ir


            # return object
            return(object)
          }
) # }}}
