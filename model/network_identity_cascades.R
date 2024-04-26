library(data.table)
library(sp)
library(foreach)
library(raster)

make_id_variables_sub <- function(id_variables){
    id_variables_sub = list()
    sub_id_variables <- function(n, id_variables_part){
        if(class(id_variables_part)=="list"){
            id_variables_sub[[n]] <<- names(id_variables_part)
            for(n1 in names(id_variables_part)){
                sub_id_variables(n1, id_variables_part[[n1]])
            }

        } 
    }
    sub_id_variables("all",id_variables)
    id_variables_sub <- id_variables_sub[names(id_variables_sub)!="all"]
    id_variables_sub <- rev(id_variables_sub)
    id_variables_sub
}

make_id_variables_full <- function(id_variables){
    id_variables_full = list()
    flatten_id_variables <- function(n, id_variables_part){
        if(class(id_variables_part[[n]])=="list"){
            for(n1 in names(id_variables_part[[n]])){
                flatten_id_variables(n1, id_variables_part[[n]])
            }
        } 
        else{ 
            id_variables_full[[n]] <<- id_variables_part[[n]]
        }
    }
    for(n in names(id_variables)) flatten_id_variables(n, id_variables)
    id_variables_full
}

initialize_data_for_cascades <- function(userFilename, networkFilename, cascadeFilename, id_variables){

    message(paste(Sys.time(), "Read user dataframe"))
    users <- readRDS(userFilename)
    setkey(users, id)
    nrow(users)
    names(users)
    head(users)
    users <<- users

    message(paste(Sys.time(), "Read edgelist dataframe"))
    edgeList <- readRDS(networkFilename)
    nrow(edgeList)
    names(edgeList)
    head(edgeList)
    Sys.time()
    edgeList <<- edgeList
    
    message(paste(Sys.time(), "Load cascade parameters"))
    cascades_info <<- readRDS(cascadeFilename)
    cascades_info$stickiness <<- 1:10/10
    
    message(paste(Sys.time(), "Fix identity variables"))
    id_variables <<- id_variables
    id_variables_sub <<- make_id_variables_sub(id_variables)
    id_variables_full <<- make_id_variables_full(id_variables)
    
    message(paste(Sys.time(), "Get Baseline Distance Distribution for Quantile Calculation"))
    N = 1000000
    set.seed(2020)
    d1 <- users[sample(users$id, N, replace=T), list(id,long,lat)]
    d2 <- users[sample(users$id, N, replace=T), list(id,long,lat)]
    if(sum(d1$id==d2$id) > 0) message("baseline distance distribution uses distances between some same users")
    d <- pointDistance(matrix(c(d1$long, d1$lat),ncol=2), 
                       matrix(c(d2$long, d2$lat),ncol=2),
                       lonlat=T)
    d_thr <- as.numeric(quantile(d, 0.25))
    rm(d1,d2)

    message(paste(Sys.time(), "Calculate identity distribution quantiles"))
    qts <- data.frame(q = seq(0,1,0.01))
    if("dist" %in% names(id_variables_full)){
        qts$dist = quantile(d, qts$q)
    } 
    for(n in unlist(id_variables_full)[unlist(id_variables_full)!='']){
        qts$var <- quantile(users[[n]], qts$q)
        names(qts)[names(qts)=="var"] <- n
    }
    qts <<- qts
    
}


load_functions_and_trial_info <- function(){
    
    get_weights <<- function(adopters, threshold_weights){
        # assign identity categories based on initial adopters
        
        # get initial adopter's identity data
        user_info <- users[adopters,]
        setkey(user_info,id)
        u1 <- rep(user_info$id, nrow(user_info))
        u2 <- rep(user_info$id, each=nrow(user_info))
        cond <- u1!=u2
        u1 <- u1[cond]
        u2 <- u2[cond]

        weights <- data.frame(w = 1)

        # calculate distances among adopters + quantile
        d <- pointDistance(matrix(c(user_info[u1,(long)], user_info[u1,(lat)]),ncol=2), 
                           matrix(c(user_info[u2,(long)], user_info[u2,(lat)]),ncol=2),
                           lonlat=T)
        if("dist" %in% names(id_variables_full)){
            weights$dist = 1 - min(qts$q[qts$dist>=median(d, na.rm=T)],na.rm=T)
        }

        # calculate differences in identity among adopters + quantile
        for(v in setdiff(names(id_variables_full),"dist")){
            qt_users <- c()
            for(n in id_variables_full[[v]]){
                qt_users <- c(qt_users,min(qts$q[qts[,n]>=median(user_info[[n]], na.rm=T)],na.rm=T))
            }
            weights$var <- max(qt_users)
            names(weights)[names(weights)=="var"] <- v
        }
        for(v in names(id_variables_sub)){
            qt_users <- c()
            for(n in id_variables_sub[[v]]){
                qt_users <- c(qt_users,weights[[n]])
            }
            weights$var <- max(qt_users)
            names(weights)[names(weights)=="var"] <- v
        }

        # choose variables above threshold_weights quantile
        wid <- weights[,names(id_variables)]
        wid_bin <- ifelse(wid>=threshold_weights,1,0)

        # if none are above threshold_weights quantile, choose the highest variable within reason
        x = threshold_weights
        while(x>threshold_weights-0.2){
            x = x - 0.05
            if(sum(rowSums(wid_bin)==0)>0)
                wid_bin[rowSums(wid_bin)==0,] <- ifelse(wid[rowSums(wid_bin)==0,]>=x,1,0)
        }
        wid_bin <- as.data.frame(wid_bin)
        wid_bin
    }
    
    find_word_identity <<- function(user_info, weights, threshold_weights){
        #assign word parameters based on median initial adopter information
        
        word_params <- list()
        XX=1

        if("dist" %in% names(id_variables)){
            long=median(user_info$long,na.rm=T)
            lat=median(user_info$lat,na.rm=T)
            users$dist_norm_w <- pointDistance(matrix(c(users$long,users$lat),ncol=2), 
                                               matrix(c(rep(long,nrow(users)), rep(lat,nrow(users))),ncol=2),
                                               lonlat=T)
        }
        for(v in setdiff(names(id_variables_full),"dist")){
            for(n in id_variables_full[[v]]){
                lab = list()
                lab[[n]] =  min(qts$q[qts[,n]>=median(user_info[[n]],na.rm=T)],na.rm=T)
                word_params[[n]] = 0
            } 
            for(n in names(lab)[lab==max(unlist(lab),na.rm=T) | unlist(lab)>threshold_weights]) word_params[[n]] = 1

            users$var <- 0
            for(n in id_variables_full[[v]]){
                users$var <- users$var + abs(users[[n]]-word_params[[n]])
            }
            names(users)[names(users)=="var"] <- paste0(v,"_norm_w")
        }
        for(v in names(id_variables_sub)){
            users$var <- 0
            for(n in id_variables_sub[[v]]){
                users$var <- users$var + users[[paste0(n,"_norm_w")]] * XX
            }
            users$var <- users$var / XX
            names(users)[names(users)=="var"] <- paste0(v,"_norm_w")
        }
        for(v in c(names(id_variables_full),names(id_variables_sub))){
            users[[paste0(v,"_norm_w")]] <- log(users[[paste0(v,"_norm_w")]] + 1)
            max_v=max(users[[paste0(v,"_norm_w")]], na.rm=T)
            min_v=min(users[[paste0(v,"_norm_w")]], na.rm=T)
            users[[paste0(v,"_norm_w")]] <- 1-(users[[paste0(v,"_norm_w")]]-min_v)/(max_v-min_v)
            word_params[[paste0("max_",v,"_norm_w")]] = max_v
            word_params[[paste0("min_",v,"_norm_w")]] = min_v
        }

        # make variable for different between user and word identity
        users$adjustment_w <<- 1
        if(rowSums(weights[,names(id_variables)])==1){
            users$adjustment_w <- 0
            for(n in names(id_variables)) users$adjustment_w <- users$adjustment_w + weights[,n]*users[[paste0(n,"_norm_w")]]
            users$adjustment_w <<- pmin(1, (users$adjustment_w-min(users$adjustment_w, na.rm=T))/(quantile(users$adjustment_w,0.95, na.rm=T) - min(users$adjustment_w, na.rm=T)))
        }
        word_params[["adjustment_w"]] <- users$adjustment_w

        word_params
    }
    
    update_edge_weights <<- function(weights, run){
        # make edge weights based on the initial adopters

        # update edge weights so they're a function of the specified parameters
        non_acc = rep(0, nrow(edgeList))
        acc = rep(0, nrow(edgeList))
        for(n in names(id_variables)){
            non_acc = non_acc + weights[,n]*edgeList[[paste0(n,"_norm")]]
            acc = acc + weights[,n]*edgeList[[paste0("acc_",n,"_norm")]] #audience accommodation (which didn't make it into the paper)
        }
        edgeList$adjustment_net <<- ifelse(rep(weights$edge,nrow(edgeList)),edgeList$weight,1) *
                                      ((1-weights$acc)*non_acc + weights$acc*acc)
        rm(non_acc, acc)

        # normalize edge weights, so the sum of each user's outdegree is 1
        edgeList <- edgeList[, names(edgeList) != "adj_total", with=F]
        edgeList <- merge(edgeList,
                          edgeList[,list(adj_total=sum(adjustment_net,na.rm=T)),by=list(to=to)],
                          by="to", all.x=T, all.y=F)
        edgeList$adjustment_net <- ifelse(edgeList$adj_total==0, 0, edgeList$adjustment_net/edgeList$adj_total)

        edgeList$adjustment <- pmin(1,edgeList$adjustment_net)
        edgeList$adjustment_original <- edgeList$adjustment
        setkey(edgeList, from, to)

        edgeList
    }
    
    update_prob <<- function(edges, exp, stickiness, thr=100){
        # update each user's probability of adoption at each iteration
        
        order <- exp$id

        f <- function(x){
            (cos(pmin(x,thr)/thr*pi)+1)/2
        }

        setkey(edges, to)    
        edges <- edges[, list(c=sum(adjustment, na.rm=T)), 
                           by = list(id = to)]
        setkey(edges, id)    
        setkey(exp, id)    
        exp[edges,p:=pmin(1,stickiness*c*adjustment_w*f(n_exposed))]

        exp[order,p]
    }

    simulate_cascade <<- function(
        artifact, adopters, num_init_adopters, 
        stickiness_vector = 0.4, p_decay_vector = 0.4, theta_vector = 100, quantile_cutoff = 0.85,
        numIterations = 100, numRuns = 1, 
        seeds = 1000, precompute_weights = FALSE, run = 1, outdir = ''
    ){

        message(paste(Sys.time(), "Initialize network"))
        adopters <- adopters[adopters %in% users$id][1:num_init_adopters]
        #if(length(adopters)==0) stop("There are no adopters in this simulation. Check that the adopters you passed in are in the users dataframe.")

        # Which identity components to include, based on initial adopters
        if(precompute_weights){
            # already precomputed in the cascades_info file
            weights = data.frame(edge = TRUE, acc = 0)
            for(v in names(id_variables)){
                weights$var = cascades_info[run,v]
                names(weights)[names(weights)=="var"] <- v
            }
        } else{
            # if you want to just pass in the adopters instead of precomputing the weights
            weights <<- cbind(data.frame(edge = TRUE, acc = 0),
                              get_weights(adopters=adopters, threshold_weights=quantile_cutoff))
        }

        #edge weights, based on initial adopters
        edgeList <- update_edge_weights(weights, run)

        #word identity, based on initial adopters
        user_info <- users[adopters,]
        word_id <- find_word_identity(user_info, weights, quantile_cutoff)

        for(which_run in 1:numRuns){ 
        for(p_decay in p_decay_vector){ 
        for(theta in theta_vector){
        for(stickiness in stickiness_vector){
            seed = seeds[which_run]

            message(paste(Sys.time(),"Run", run, artifact, "Adopters", length(adopters)))

            # Make Filename
            name = c(artifact, "-")
            if(weights$edge) name = c(name, "weight-") #tie strength included
            for(n in names(id_variables)) name = c(name, n, round(weights[,n],2), "-") #which id variables included
            name = c(name, # other parameters
                     "p_decay", round(p_decay,2), "-",
                     "thr", theta, "-",
                     "S", round(stickiness,2))
            name = paste0(name, collapse="")
            message(name)

            #Initialize simulation
            set.seed(seed)
            iteration = 1; num_adopt = 0
            exposed <- data.table(id=adopters, 
                                  by="initial", p=1, first_exposed=0, n_exposed=0,
                                  first_adopt=1, n_adopt=0)
            setkey(exposed,id)
            exposed$adjustment_w <- users[exposed$id,(adjustment_w)]
            edgeList$adjustment <- edgeList$adjustment_original
            ts_sim <- data.frame()

            #Start Runs
            message(paste(Sys.time(),"start iterations"))
            while(iteration < numIterations | iteration %% 10 !=0 | 
                  sum(exposed$n_adopt>0)>=1.01*num_adopt){ #keep going until under 1% increase in adoption over 10 iterations
                if(iteration %% 10 == 0) num_adopt <- sum(exposed$n_adopt>0)
                if(iteration %% 50 == 0) message(paste(Sys.time(),"decide who uses the new word"))

                #each agent decides whether to use the word
                use <- runif(nrow(exposed),0,1) < exposed$p
                use[is.na(use)] <- FALSE
                if(sum(use)>0){
                    if(iteration %% 50 == 0) message(paste(Sys.time(),"update adopt list"))

                    #update adpoption information (n_adopt, first_adopt) for agents who adopted the artifact this iteration
                    exposed[use, n_adopt := n_adopt + 1]
                    exposed[use & is.na(first_adopt), first_adopt := iteration] 

                    #update adoption timeseries dataframe for agents who adopted the artifact this iteration
                    adopted <- exposed$n_adopt > 0
                    ts_sim <- rbind(ts_sim, edgeList[exposed$id[adopted],]
                                                    [to %in% exposed$id[use], 
                                                     list(neigh=paste(from,adjustment,sep=":",collapse=";"), 
                                                          iteration=iteration), by = list(id = to)])

                    #find and refresh exposed edges
                    if(iteration %% 50 == 0) message(paste(Sys.time(),"update exposure list"))
                    exp_id <- sort(unique(edgeList[exposed$id[use],(to)]))
                    edgeList[exp_id, adjustment := adjustment_original] 
                    edges <- edgeList[exposed$id[use],]

                    ##for agents exposed previously: update n_exposed
                    prev <- edges$to %in% exposed$id
                    if(sum(prev)>0){
                        e <- edges[prev, list(nexp=length(to)), by=list(id=to)]
                        exposed[e$id, n_exposed := n_exposed + e$nexp]
                    } 

                    ##for agents' first exposure: append rows to exposed corresponding to these agents
                    if(sum(!prev)>0) {
                        e <- edges[!prev, list(by=paste(from,collapse=","),
                                           p=0, first_exposed=iteration,
                                           n_exposed=length(to),
                                           first_adopt=NA, n_adopt=0, adjustment_w=0), 
                                   by=list(id=to)]
                        e$adjustment_w <- users[e$id,(adjustment_w)]
                        exposed <- rbindlist(list(exposed, e))
                    }
                    setkey(exposed, id)

                    #update usage probability of exposed agents
                    if(iteration %% 50 == 0) message(paste(Sys.time(),"update probability"))
                    if(sum(exposed$n_adopt>0) > length(adopters)) exposed[, p := (1-p_decay)*p] # decay non-exposed edges
                    exposed[exp_id, 
                            p := update_prob(edges, # update probability for exposed users edges
                                             exposed[exp_id,list(id,p,adjustment_w,n_exposed)], 
                                             stickiness=stickiness,
                                             thr=theta)]
                    if(iteration %% 50 == 0) message(paste(Sys.time(),"iteration done"))
               }

                if(nrow(exposed)>length(adopters)) iteration = iteration + 1
                if(iteration %% 50 == 1) message(paste(Sys.time(), iteration, 
                                                       nrow(exposed), sum(exposed$n_adopt>0)))
            }

            # Save Outputs

            message(paste(Sys.time(),"aggregate by geography"))
            exposed <- merge(users, exposed, by="id", all.x=T, all.y=F)
            setkey(exposed,id)
            r <- exposed[,list(num_users = length(id),
                                 exposed.0.=quantile(first_exposed,0, na.rm=T),
                                 exposed.25.=quantile(first_exposed,0.25, na.rm=T),
                                 adopt.0.=quantile(first_adopt,0, na.rm=T),
                                 adopt.25.=quantile(first_adopt,0.25, na.rm=T),
                                 n_exposed.total=sum(n_exposed,na.rm=T),
                                 n_exposed.median=median(n_exposed,na.rm=T),
                                 n_exposed.atall=sum(n_exposed>0,na.rm=T),
                                 n_adopt.total=sum(n_adopt,na.rm=T),
                                 n_adopt.median=as.double(median(n_adopt,na.rm=T)),
                                 n_adopt.atall=sum(n_adopt>0,na.rm=T)
                              ),
                           by=list(GEOID=geo_fips, state=state, division=division, region=region)] 
            r$seed <- seed
            r$artifact <- artifact

            message(paste(Sys.time(),"save tables"))
            saveRDS(exposed, paste0(outdir, "exposed-",
                                    name, "-", run, "-", seed, ".RDS"))  
            saveRDS(ts_sim, paste0(outdir, "ts-", 
                                   name, "-", seed, ".RDS"))
            saveRDS(r, paste0(outdir, "geo-", 
                              name, "-", seed, ".RDS"))
        }
        }
        }
        }
    }
}