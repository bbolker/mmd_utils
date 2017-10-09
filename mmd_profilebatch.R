load("allfits.RData")
profList <- list()
names(profList) <- names(allfits)
for (i in seq_along(allfits)) {
    cat(names(allfits)[i],"\n")
    profList[[i]] <- profile(allfits[[i]])
    save("profList",file="allprofs.RData")
}

