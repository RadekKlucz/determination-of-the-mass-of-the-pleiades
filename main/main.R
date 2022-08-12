# Load mass of a star cluster data 
path <- getwd()
setwd(path)
mass <- read.table('./data/pleiades.txt', header = TRUE)

# Analyzing data of mass star cluster 
summary(mass)

# Defining the function using in the script 
compute_degrees <- function(data, prefix) { 
  if(prefix == 'RAJ'){ 
    right_ascension =  data.frame(data$RAJ2000_1 + data$RAJ2000_2 / 60 + data$RAJ2000_3 / 360)
    return(right_ascension)    
    }
  else { 
    declination <- data.frame(data$DEJ2000_1 + data$DEJ2000_2 / 60 + data$DEJ2000_3 / 360)
    }
}

generall_function <- function(data, prefix) {
  distance <- data.frame(1 / (data$plx * 10 ^ -3))
  pc_dec <- data.frame(distance * data$pmDE * 4.84814e-6 * 10 ^ -3)
  pc_ra <- data.frame(distance * data$pmRA * 4.84814e-6 * 10 ^ -3)
  if(prefix == 'plx'){
    return(distance)    
  }
  else if(prefix == 'pmDE') { 
    names(pc_dec)[1] <- 'dec' 
    return(pc_dec)
  }
  else { 
    names(pc_ra)[1] <- 'ra' 
    return(pc_ra)
  }
}

# In hours
right_ascension <- compute_degrees(mass, 'RAJ')
names(right_ascension)[1] <- 'ra'

# In degrees
declination <- compute_degrees(mass, 'DEJ')
names(declination)[1] <- 'dec'

# Defining variable g constant of gravity
g = 6.6743 * 10  ^ -11

# Calculate distance, right ascension and declination using generall function 
distance <- generall_function(mass, 'plx')
names(distance)[1] <- 'dist'
pc_declination <- generall_function(mass, 'pmDE')
pc_right_ascension <- generall_function(mass, 'pmRA')

head(distance)
head(pc_declination)
head(pc_right_ascension)

# Calculate average values in pc/yr
average_distance <- mean(distance$dist)
error_distance <- sd(distance$dist / length(distance$dist))
average_pc_ra <- mean(pc_right_ascension$ra)
average_pc_dec <- mean(pc_declination$dec)
ra_center <- mean(right_ascension$ra)
dec_center <- mean(declination$dec)

pm_ra_cor <- pc_right_ascension - average_pc_ra
head(pm_ra_cor)
pm_dec_cor <- pc_declination - average_pc_dec
head(pm_dec_cor)

# Calculate in meters per second
pm_ra_cor_mps <- pm_ra_cor * 977813106
pm_dec_cor_mps <- pm_dec_cor * 977813106
head(pm_ra_cor_mps)
head(pm_dec_cor_mps)

# Calculate velocity in meters per second
velocity_sqr <- pm_ra_cor_mps ^ 2 + pm_dec_cor_mps ^ 2
head(velocity_sqr)

# Calculate averge distance
average_distance_m <- average_distance * 3.08567758128e+16

# Calculate average velocity in meters per second ^ 2
average_vel_sqr <- mean(velocity_sqr$ra) * (3 / 2)

# Calculate angel size 
angle_size <- max(declination$dec) - min(declination$dec)
size <- average_distance * angle_size * 0.0174533

# Calculate angel size in meters
size_meter <- size * 3.086e+16

# Calculate mass a star cluster
mass_cluster <- (size_meter / 2 * average_vel_sqr) / g

# Calculate your own star cluster movements
movements_ra <- mean(pm_ra_cor$ra)

movements_dec <- mean(pm_dec_cor$dec)

# Calculate distance stars from a star cluster
deg_to_rad <- 0.0174533
cluster <- data.frame(sqrt(distance$dist ^ 2 + average_distance ^ 2 - 2 * distance$dist * 
                  average_distance * (sin(declination$dec * deg_to_rad) * sin(dec_center * 
                  deg_to_rad) * cos((right_ascension$ra - ra_center) * deg_to_rad) +
                  cos(declination$dec * deg_to_rad) * cos(dec_center * deg_to_rad))))
names(cluster)[1] <- 'cluster-r'

# Create data frame 

data <- data.frame(ra=right_ascension, dec=declination, dec_center=dec_center, 
                   average_pm_ra=average_pc_ra, average_pm_dec=average_pc_dec,
                   pm_ra=pm_ra_cor, pm_dec=pm_dec_cor, pm_rec_cor_i=pm_ra_cor_mps, 
                   pm_dec_cor_i=pm_dec_cor_mps, velocity=velocity_sqr, 
                   dist_center=average_distance, cluster_r=cluster)
head(data)

# Sort values
data[order(data$cluster.r),]
sorted_velocity <- data$ra.3

# Split data for 60 parts
nr <- length(sorted_velocity)
parts <- split(sorted_velocity, rep(1:ceiling(nr / 19), each=19, length.out=nr))

# Create list of quantile and calculate average value
list_of_quantile = list()
for(i in parts) { 
  list_of_quantile <- append(list_of_quantile, quantile(i, probs = c(.90)))
}

average_list <- mean(list_of_quantile$`90%`)

# Calculate ending mass
mass_cluster_end <- (size_meter / 2 * average_list) / g

# Print all needed values
print(paste0('Average distance: ', average_distance))
print(paste0('Linear velocities of stars in right ascension and declination: ', 
             average_pc_ra,' ',  average_pc_dec, ' pc'))
print(paste0('Separated own motions of individual stars, average values: ', 
             movements_ra, ' ',  movements_dec, ' pc/yr'))
print(paste0('Mean square velocity value: ', average_vel_sqr, ' m/s'))
print(paste0('Cluster angular size: ', size))
print(paste0('Mass: ', mass_cluster, ' kg'))
