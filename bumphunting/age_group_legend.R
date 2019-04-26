library('RColorBrewer')

pdf('pdf/age_groups.pdf', useDingbats = FALSE)
plot.new()

labs <- c(paste(c('Glia:', 'Neuron:'), rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2)), 'Homogenate Prenatal')
cols <- c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50')
legend('top', labs[c(1, 3, 5, 7, 9, 2, 4, 6, 8)], col = cols[c(1, 3, 5, 7, 9, 2, 4, 6, 8)], lwd = 3, bty = 'n', ncol = 2)
dev.off()

pdf('cell_types.pdf', useDingbats = FALSE)
plot.new()
legend(x = 'top', c("Glia", "Neuron"), col = brewer.pal(name="Dark2", n=2), lwd = 3, bty = 'n', ncol = 2)
dev.off()

pdf('./bumphunting/pdf/age_groups_1column.pdf', useDingbats = FALSE)
plot.new()
labs <- c(paste(c('Glia:', 'Neuron:'), rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2)), 'Homogenate Prenatal')
cols <- c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50')
legend('top', labs[c(2, 4, 6, 8, 1, 3, 5, 7, 9)], fill = cols[c(2, 4, 6, 8, 1, 3, 5, 7, 9)], bty = 'n', ncol = 1)
dev.off()

pdf('./bumphunting/pdf/age_groups_1column_line.pdf', useDingbats = FALSE)
plot.new()
labs <- c(paste(c('Glia:', 'Neuron:'), rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2)), 'Homogenate Prenatal')
cols <- c(brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)], 'grey50')
legend('top', labs[c(2, 4, 6, 8, 1, 3, 5, 7, 9)], col = cols[c(2, 4, 6, 8, 1, 3, 5, 7, 9)], lwd = 3, bty = 'n', ncol = 1)
dev.off()


pdf('./bumphunting/pdf/age_groups_1column_noPrenatal.pdf', useDingbats = FALSE)
plot.new()
labs <- paste(c('Glia:', 'Neuron:'), rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2))
cols <- brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)]
legend('top', labs[c(1, 3, 5, 7, 2, 4, 6, 8)], fill = cols[c(1, 3, 5, 7, 2, 4, 6, 8)], bty = 'n', ncol = 1)
dev.off()

pdf('./bumphunting/pdf/age_groups_4column_noPrenatal.pdf', useDingbats = FALSE)
plot.new()
labs <- paste(c('Glia:', 'Neuron:'), rep(c('Infant', 'Child', 'Teen', 'Adult'), each = 2))
cols <- brewer.pal(8, "Paired")[c(5:6, 7:8, 3:4, 1:2)]
legend('top', labs, fill = cols, bty = 'n', ncol = 4)
dev.off()


