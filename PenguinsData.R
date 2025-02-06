# Call the palmer penguins data
dat.full <- palmerpenguins::penguins
#dat

# Number of missing values
na.num <- sum(is.na(dat.full))
#na.num

# Number of missing values in variable subset
na.sub.num <- sum(is.na(dat.full[,c(1,3,4,5,6)]))

# Identify which observations contain missing values
miss.ind <- which(apply(dat.full[,c(1,3,4,5,6)],1,function(x)any(is.na(x))))
#miss.ind

# Investigate the two observations with missing values
miss.vals <- dat.full[miss.ind,c(1,3,4,5,6)]
# miss.vals

# Palmer Penguins data Variable subset with missing observations removed
dat <- dat.full[-miss.ind,c(1,3,4,5,6)]

# Colour Palette for the penguins
Class.int <- as.numeric(unlist(dat[,1]))
col.vec <- adjustcolor(RColorBrewer::brewer.pal(n=3, name='Dark2'),alpha.f = 0.5)
cols <- col.vec[Class.int]
