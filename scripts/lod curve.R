# Load required library
library(ggplot2)

# Define LOD score data
locus_position <- c(0, 25, 40, 55, 65, 75, 90, 100)  # Positions of markers
lod_score <- c(1.5, 4.0, 5, 3.0, 1.0, 2.5, 3.0, 1.5)  # Corresponding LOD scores

# Create a data frame
df <- data.frame(locus_position, lod_score)

# Define the significance threshold
lod_threshold <- 3.5

# Plot the LOD score curve with smooth interpolation
ggplot(df, aes(x = locus_position, y = lod_score)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.4, aes(color = "LOD Curve"), size = 2) +  # Smooth curve
  geom_hline(yintercept = lod_threshold, linetype = "dashed", color = "red",size=2) +  # Threshold line with legend+  # Threshold line
  geom_point(x=40,y=5,color="blue",size=4)+
  #annotate("text", x = 50, y = 5.5, label = "Maximum likelihood QTL between loci G and H", size = 10) +
  geom_segment(aes(x = 55, xend = 42, y = 5, yend = 5), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               size = 1) +  # Arrow pointing to peak +
  scale_color_manual(values = c("LOD Curve" = "black")) +  # Customize colors
  theme_minimal() +
  labs(x = "Locus Position", y = "LOD Score",size= 5) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), 
                     labels = c("0", "2", "4", "8", "10", "12", "14", "16", "18", "19", "X")) +  # Matching labels
  scale_y_continuous(breaks = seq(0, 5.5, by = 1)) +
  theme_minimal() +
  theme(
    legend.position = "none", 
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    #axis.line = element_line(color = "black", size = 1),  # Add x and y axes
    panel.grid.major = element_blank(),  # Remove major grid lines for clarity
    panel.grid.minor = element_blank(),# Remove minor grid lines
    axis.title.x = element_text(size = 24),  # Increase x-axis label font size
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20),   # Resize x-axis values
    axis.text.y = element_text(size = 20)
  )

print("DAIR3 Project")
