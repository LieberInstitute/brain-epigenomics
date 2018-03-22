library('RColorBrewer')

pdf('pdf/age_groups.pdf', useDingbats = FALSE)
plot.new()

labs <- c(paste(c('Glia:', 'Neuron:'), rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2)), 'Homogenate Prenatal')
cols <- c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50')
legend('top', labs[c(1, 3, 5, 7, 9, 2, 4, 6, 8)], col = cols[c(1, 3, 5, 7, 9, 2, 4, 6, 8)], lwd = 3, bty = 'n', ncol = 2)
dev.off()
