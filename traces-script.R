
# This script takes fsa files from sequencing
# It aligns the ladder to the sample trace
# converts from scan size to bp
# and exports the sizing data


# install packages if not already installed
# install.packages("BiocManager") 
# BiocManager::install("sangerseqR")
# install.packages("pracma") 

# load libraries
library(sangerseqR)
library(pracma)

# set working directory. use " " Use / Not \
setwd("")

fsa <- read.abif("A07_ME20250911_nuclease2_A07_2025-09-11_1.fsa")

# Access a dye channel (blue trace is in data 1)
trace <- fsa@data$DATA.1

# view the raw fluorescence signal
plot(trace, type = "l", col = "blue", main = "Electropherogram Trace", xlab = "Scan", ylab = "Fluorescence")

# ladder trace
ladder_trace <- fsa@data$DATA.105

# view liz 500 ladder
plot(fsa@data$DATA.105, type = "l", col = "red", main = "Ladder Trace (DATA.105)")


###################################################################################################

#find peaks
peaks <- findpeaks(ladder_trace, threshold = 200, minpeakdistance = 20) # can change threshold/minpeakdistance

#liz500 ladder - 16 peaks
known_sizes_bp <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)

scan_positions <- peaks[,2]  # peak locations

#################################################################################
# if A MISMATCH WAS FOUND BETWEEN NUMBER OF
# PEAKS AND NUMBER OF KNOWN LADDER PEAKS
# CAN CHECK:
length(known_sizes_bp)
length(scan_positions)

# CHECK ITS CALLED RIGHT PEAKS
plot(ladder_trace, type = "l", main = "Ladder Trace")
points(scan_positions, ladder_trace[scan_positions], col = "red", pch = 19)


#if mismatch, need to correct it, otherwise go to next section

# IF THERE ARE EXTRA PEAKS - IDENTIFY THEIR SIZES:
sort(scan_positions)

# Manually choose the correct peak indices
# exclude the wrong ones
# correct_scan_positions <- scan_positions[2:17]
# or
# exclude the ones in bracket
correct_scan_positions <- scan_positions[!(scan_positions %in% c(1248, 1360, 2878))]

#Now both vectors have same length: check
length(correct_scan_positions)  # Should be 16(liz500) 
length(known_sizes_bp)          # Should be 16(liz500) 

#FIXED LADDER
plot(ladder_trace, type = "l", main = "Ladder Trace")
points(correct_scan_positions, ladder_trace[correct_scan_positions], col = "red", pch = 19)

#CORRECT THE SIZES if need to
#known_sizes_bp <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)  # example ladder sizes

#end of correcting ladder section
#############################################################################


# Plot ladder trace with peaks marked
plot(ladder_trace, type = "l", col = "orange", main = "FAM & Ladder Traces with Peaks",
     xlab = "Scan Position", ylab = "Fluorescence Intensity", ylim = range(c(trace, ladder_trace)))
points(correct_scan_positions, ladder_trace[correct_scan_positions], col = "red", pch = 19)
# Add fam trace on same plot but scaled and in orange
lines(trace, col = "blue")
legend("topright", legend = c("Ladder", "Ladder Peaks", "FAM Trace"),
       col = c("orange", "red", "blue"), lty = c(1, NA, 1), pch = c(NA, 19, NA))

###########################################################################

# order the peaks
correct_scan_positions_sorted <- sort(correct_scan_positions)
  
df_peaks <- data.frame(scan = correct_scan_positions_sorted, size = known_sizes_bp)
  
############################
# linear
calibration <- lm(size ~ scan, data = df_peaks)

# scan number to bp conversion
scan_to_bp <- function(scan_pos) {
   coef(calibration)[1] + coef(calibration)[2] * scan_pos
  }
  
# Apply to entire traces
all_ladder_scans <- seq_along(ladder_trace)
all_ladder_bp <- scan_to_bp(all_ladder_scans)
  
all_fam_scans <- seq_along(trace)
all_fam_bp <- scan_to_bp(all_fam_scans)
  
  
# Plot an empty plot with correct limits (xlim and ylim)
plot(NULL, xlim = c(0, 80), ylim = c(0, 5000),
      main = "FAM Trace",
       xlab = "Base Pair (bp)", ylab = "Fluorescence Intensity")
  
# Add fam trace in blue
  lines(all_fam_bp, trace, col = "blue")
  
  
 
  
#######################
# Create data frame for ladder nuclease
ladder_df <- data.frame(
    Scan = all_ladder_scans,
    BasePair = all_ladder_bp,
    Fluorescence = ladder_trace
  )
  
# Create data frame for FAM nuclease
fam_df <- data.frame(
    Scan = all_fam_scans,
    BasePair = all_fam_bp,
    Fluorescence = trace
  )
  
  
# Filter FAM trace for 0–200 bp - you can change range
fam_df_subset <- subset(fam_df, BasePair >= 0 & BasePair <= 200) 
 
# Export to CSV files NUCLEASE
write.csv(ladder_df, "ladder_trace_with_bp.csv", row.names = FALSE)
write.csv(fam_df_subset, "fam_trace_with_bp.csv", row.names = FALSE)
  

######CAN OVERLAY TRACES IN GRAPHPAD#####