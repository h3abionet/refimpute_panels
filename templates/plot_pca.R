#!/usr/bin/env Rscript
# Title     : TODO
# Objective : TODO
# Created by: mamana
# Created on: 2019-06-21

library(ggplot2)
library(reshape2)

shapes <- c(0:25)
evecDat = read.table("${input_evec}", col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "Pop", "Group", "Supergroup"))

if (nlevels(evecDat\$Group) > 1) {
    output_pdf <- "${output_pdf}_pc1_pc2.tiff"
    p <- ggplot(data=evecDat, aes(x=PC1, y=PC2, colour=Supergroup, shape=Supergroup)) +
    geom_point(size=2) +
    scale_shape_manual(values=rep(shapes, times=4))
    p + theme( panel.background= element_rect(fill="#E3E3EE"))
    ggsave(output_pdf, units='in', dpi=120)

    # scale_shape_identity(breaks = unique(evecDat\$Pop), guide = 'legend')
    output_pdf <- "${output_pdf}_pc1_pc3.tiff"
    p <- ggplot(data=evecDat, aes(x=PC1, y=PC3, colour=Supergroup, shape=Supergroup)) +
        geom_point(size=2) +
        scale_shape_manual(values=rep(shapes, times=4))
    p + theme( panel.background= element_rect(fill="#E3E3EE"))
    ggsave(output_pdf, units='in', dpi=120)
} else {
    shapes <- c(66,68,71,18,0,1,2,3,4,5,6,7,8,35,64,9,36,10,11,65,66,12,13,67)
    output_pdf <- "${output_pdf}_pc1_pc2.tiff"
    p <- ggplot(data=evecDat, aes(x=PC1, y=PC2, colour=Supergroup, shape=Supergroup)) +
    geom_point(size=2) +
    scale_shape_manual(values=rep(shapes, times=4))
    p + theme( panel.background= element_rect(fill="#E3E3EE"))
    ggsave(output_pdf, units='in', dpi=120)

    output_pdf <- "${output_pdf}_pc1_pc3.tiff"
    p <- ggplot(data=evecDat, aes(x=PC1, y=PC3, colour=Supergroup, shape=Supergroup)) +
    geom_point(size=2) +
    scale_shape_manual(values=rep(shapes, times=4))
    p + theme( panel.background= element_rect(fill="#E3E3EE"))
    ggsave(output_pdf, units='in', dpi=120)
}
