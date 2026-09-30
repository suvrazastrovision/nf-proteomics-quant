#!/usr/bin/env Rscript
library(plotly)
library(glue)
library(tidyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--precursors', type='character', nargs='+')
parser$add_argument('--non-norm-precursors', type='character', nargs='+', default=c())
parser$add_argument('--inputfn', type='character')
parser$add_argument('--conflvl', type='double')
opt = parser$parse_args()

precursortables = opt$precursors
nn_precursors = opt$non_norm_precursors
inputfnpath = opt$inputfn
enzyme = opt$enzyme

# ModSeq, charge, file combination is "a single precursor"
seqcol = 'Modified.Sequence'
chargecol = 'Precursor.Charge'
filenamecol = 'Run'
miscleavcol = 'missed_cleavages'


# For boxplots:
rtcol = 'RT'
precerrcol = 'Ms1.Apex.Mz.Delta'
precquantcol = 'Ms1.Area'
peakwidthcol = 'FWHM'
ms2quantcol = 'Precursor.Quantity'
samplecol = 'sample'

amount_ms2 = read.table("concat_filescans", sep="\t", header=F)
colnames(amount_ms2) = c('file', 'nr_scans')

boxplot_stats = function(data, col) {
  summary_stats = data %>%
  summarise(lower = quantile({{col}}, 0.25, na.rm=T),
            middle = median({{col}}, na.rm=T),
            upper = quantile({{col}}, 0.75, na.rm=T),
            minval = min({{col}}, na.rm=T),
            maxval = max({{col}}, na.rm=T),
    )
  summary_stats$iqr = summary_stats$upper - summary_stats$lower
  summary_stats$whisk_min = pmax(summary_stats$minval, summary_stats$middle - summary_stats$iqr * 1.5)
  summary_stats$whisk_max = pmin(summary_stats$maxval, summary_stats$middle + summary_stats$iqr * 1.5)
  return(summary_stats)
}
colmap = list(
  sample=c(samplecol),
  file=c(filenamecol)
  )
## List of plot names and their axis label:
ptypes = list(
  retentiontime=c(rtcol, 'time(min)'),
  precerror=c(precerrcol, 'Precursor error (ppm)'),
  precquant=c(precquantcol, 'Quantity'),
  peakwidth=c(peakwidthcol, 'FWHM')
  # Is there no # of fragments data in DIANN?
)
  
inputfn = read.table(inputfnpath, header=T, sep="\t", comment.char = "", quote = "")
inputfn$Run = tools::file_path_sans_ext(basename(inputfn$file_path))
for (grouper in names(colmap)) {
  chunk = 0
  xcol = colmap[[grouper]][1]
  accu_miscleav = data.frame(file=character(), precursorcount=character(), sample=character(), nr_scans=character(), missed_cleavage=character(), nrprec=character(), text=character(), percent=character())
  accu_nrprec = data.frame(file=character(), precursorcount=character(), sample=character(), nr_scans=character())
  for (precursortable in precursortables) {
    precs_labeled = read.table(precursortable, header=T, sep="\t", comment.char = "", quote = "")
    nr_verts = length(unique(precs_labeled$Run))
    
    # 2200 is fine if pixels?, make that cutoff then, 100 samples (200 * 100 + 200) =
    # 2400 is not
    vert_height = 20 * nr_verts + 200
    
    samplefnmap = unique(data.frame(file=precs_labeled[[filenamecol]], sample=precs_labeled[[samplecol]]))
    precursors = aggregate(precs_labeled[c(seqcol)], by=precs_labeled[xcol], length)
    names(precursors) = c(grouper, 'precursorcount')
  
    miscleav = aggregate(precs_labeled[c(seqcol)], by=precs_labeled[c(miscleavcol, xcol)], length)
    names(miscleav) = c('missed_cleavage', grouper, 'nrprec')
    miscleav$text = glue('{miscleav$nrprec} precursors')
    # Subset so we can use this when plotting
    miscleav = subset(miscleav, missed_cleavage %in% c(0,1,2))
    if (grouper == 'file') {
      fake_3_mc = data.frame(missed_cleavage=3, file=unique(miscleav$file), nrprec=0, text=unique(miscleav$file))
      miscleav_plot = rbind(miscleav, fake_3_mc)
      precursors = merge(precursors, samplefnmap, by=grouper, all=T)
      precursors = merge(precursors, amount_ms2, by=grouper)
    } else {
      miscleav_plot = miscleav
    }
  
    miscleav = merge(precursors, miscleav, by=grouper)
    miscleav_plot = merge(precursors, miscleav_plot, by=grouper)
    miscleav$percent = miscleav$nrprec / miscleav$precursorcount * 100
    miscleav_plot$percent = miscleav_plot$nrprec / miscleav_plot$precursorcount * 100
    if (grouper == 'file') {
      mc_text_y = 50
      mc_text_size = 3
      accu_miscleav = rbind(accu_miscleav, miscleav)
      accu_nrprec = rbind(accu_nrprec, precursors)
    } else {
      miscleav_plot$percent = miscleav_plot$nrprec / miscleav_plot$precursorcount * 100
      mc_text_y = max(miscleav_plot$percent) * 2/6
      mc_text_size = 4
    }
    
    # Barplot nr of precursors
    ggp = ggplot(precursors) +
      geom_bar(aes(x=.data[[grouper]], y=precursorcount), stat='identity', position='dodge') + 
      coord_flip() + 
      ylab('# precursors') + 
      scale_x_discrete(labels=precursors[['sample']]) +
      theme_bw() + 
      theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
            axis.text.x=element_text(size=10),
            legend.text=element_text(size=10), legend.title=element_blank())
      if (grouper == 'file') {
        # plotly doesnt support hjust, so calculate position to be in middle
        ggp = ggp +
          geom_text(aes(x=.data[[grouper]], y=precursorcount / 2, label=.data[[grouper]]), size=3, colour="white")
      } else {
        ggp = ggp + theme(axis.text.y=element_text(size=10, angle=90)) +
          geom_text(aes(x=.data[[grouper]], y=precursorcount / 2, label=precursorcount), size=10, colour="white")
      }
  
    p = ggplotly(ggp, width=600, height=vert_height) %>%
            layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    # Work around since plotly does not honor above legend.title=element_blank call
    p$x$layout$legend$title$text = ''
    htmlwidgets::saveWidget(p, glue('precursorplothtml/{grouper}__{chunk}__amount_precursors.html'), selfcontained=F)
  
    ## Missed cleavages plot
    miscleav_plot$missed_cleavage = as.factor(miscleav_plot$missed_cleavage)
    mcplot = ggplot(miscleav_plot) +
        geom_bar(aes(x=.data[[grouper]], y=percent, fill=missed_cleavage, group=missed_cleavage), position='dodge', stat='identity') +
        # 0.9 is the default dodge (90% of 1, 1 used bc all same value) but when not spec -> no dodge at all?
        geom_text(position=position_dodge(width=0.9), aes(x=.data[[grouper]], y=mc_text_y, group=missed_cleavage, label=text), colour="black", size=mc_text_size, inherit.aes=T) +
        ylim(c(0, 100)) + ylab('% of precursors') +
        scale_x_discrete(labels=miscleav_plot[['sample']]) +
        theme_bw() +
        theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
              legend.position="top", legend.text=element_text(size=10), legend.title=element_blank()) +
        coord_flip() 
  
  
    p = ggplotly(mcplot, width=600, height=vert_height) %>%
            layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    p$x$layout$legend$title$text = ''
    htmlwidgets::saveWidget(p, glue('precursorplothtml/{grouper}__{chunk}__missed_cleavages.html'), selfcontained=F)
  
    for (ptype in names(ptypes)) {
      if (ptypes[[ptype]][1] %in% colnames(precs_labeled)) {
        summary_stats <- precs_labeled[c(filenamecol, ptypes[[ptype]][1])] %>%
          group_by(.data[[filenamecol]]) %>%
          boxplot_stats(.data[[ptypes[[ptype]][1]]])
    
        summary_stats$Run = as.character(summary_stats$Run)
        sumstats_samples = merge(summary_stats, inputfn[,c('Run', samplecol)], by='Run')
    
        # Plot using geom_crossbar, geom_errorbar
        ggp = ggplot(sumstats_samples, aes(x=Run, y=middle)) +
          geom_errorbar(aes(ymin = whisk_min, ymax = whisk_max), width=0) +
          geom_crossbar(aes(ymin = lower, ymax=upper), fill='white', linewidth=0.15) +
          scale_x_discrete(labels=sumstats_samples$sample) +
          coord_flip() + theme_bw() +
          ylab(ptypes[[ptype]][2]) + 
          theme(axis.title=element_text(size=15), axis.title.y=element_blank(),
            axis.text.x=element_text(size=10), axis.text.y=element_text(angle=90)) +
          # plotly doesnt support hjust, so calculate position to be in middle
          geom_text(aes(x=Run, y=(get('whisk_min') + get('whisk_max')) / 2, label=Run), position=position_nudge(x=0.1), size=3, colour="black")
        if (ptype == 'precerror') { ggp = ggp + geom_hline(yintercept=0, size=1, color='green3', linetype='dashed') }
        p = ggplotly(ggp, width=600, height=vert_height)
        htmlwidgets::saveWidget(p, glue('precursorplothtml/{ptype}__{chunk}.html'), selfcontained=F)
      }
    }
    chunk = chunk + 1
  }
}
write.table(accu_miscleav, glue('miscleav_qc.txt'), row.names=F, quote=F, sep='\t')
write.table(accu_nrprec, glue('counttable_qc.txt'), row.names=F, quote=F, sep='\t')

  
  
  
#### Precursor table also contains protein/gene group data, so do all in single script
# Summary tables of overlap
featmap = list(proteins=c('Protein.Group', 'PG.Q.Value', 'PG.MaxLFQ'), genes=c('Genes', 'GG.Q.Value', 'Genes.MaxLFQ'))
for (feattype in names(featmap)) {
  chunk = 0
  accumulated_feats = data.frame(feat=character())
  accumulated_feats_q = data.frame(feat=character())
  for (precursortable in precursortables) {
    precs_labeled = read.table(precursortable, header=T, sep="\t", comment.char = "", quote = "")
    featcol = featmap[[feattype]][1]
    featfiltcol = featmap[[feattype]][2]
    quantcol = featmap[[feattype]][3]
    prec_featcols = precs_labeled[, c('Run', samplecol, featmap[[feattype]]) ]
    precfeats_filt = prec_featcols[prec_featcols[featfiltcol] < opt$conflvl, ]
    # Keep only one precursor line per feature per run
    feats = precfeats_filt %>%
      distinct(Run, .data[[featcol]], .keep_all=T)
  
    # Remove NA/0 for both precursor and feats table:
    precfeats_filt = precfeats_filt[precfeats_filt[[quantcol]] >= 0,]
    feats_filtered = feats[feats[[quantcol]] >= 0,]
  
    # Summary overlap table data 
    feats_nrfns = aggregate(get(filenamecol)~get(featcol), feats, length)
    feats_nrfns_q = aggregate(get(filenamecol)~get(featcol), feats_filtered, length)
    colnames(feats_nrfns) = c('feat', 'nrfns')
    colnames(feats_nrfns_q) = c('feat', 'nrfns')
    accumulated_feats = merge(accumulated_feats, feats_nrfns, by='feat', all.x=T, all.y=T)
    accumulated_feats_q = merge(accumulated_feats_q, feats_nrfns_q, by='feat', all.x=T, all.y=T)

    # Missing features plot
    feats_missing = feats[feats[[quantcol]] <= 0,]
    if (nrow(feats_missing)) {
      missing_per_run = aggregate(get(quantcol)~Run, feats_missing, length)
      colnames(missing_per_run)[2] = 'nr_missing'
      missing_samples = merge(missing_per_run, inputfn[,c('Run', samplecol)], by='Run')
      ggp = ggplot(missing_per_run) +
        geom_bar(aes(x=Run, y=nr_missing), stat='identity', position='dodge') + 
        coord_flip() + 
        ylab('# missing features') + 
        scale_x_discrete(labels=missing_samples$sample) +
        theme_bw() + 
        theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
              axis.text.x=element_text(size=10), axis.text.y=element_text(angle=90),
              legend.text=element_text(size=10), legend.title=element_blank()) +
        # plotly doesnt support hjust, so calculate position to be in middle
        geom_text(aes(x=Run, y=nr_missing / 2, label=Run), size=3, colour="white")
      p = ggplotly(ggp, width=600, height=vert_height) %>%
              layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
      # Work around since plotly does not honor above legend.title=element_blank call
      p$x$layout$legend$title$text = ''
      htmlwidgets::saveWidget(p, glue('{feattype}plothtml/missing_feats__{chunk}.html'), selfcontained=F)
    }
   
  
    # Total features plot
    feats_per_run = aggregate(get(quantcol)~Run, feats_filtered, length)
    colnames(feats_per_run)[2] = 'nr_feats'
    feats_per_run_samples = merge(feats_per_run, inputfn[,c('Run', samplecol)], by='Run')
    ggp = ggplot(feats_per_run) +
      geom_bar(aes(x=Run, y=nr_feats), stat='identity', position='dodge') + 
      coord_flip() + 
      ylab('# identified') + 
      scale_x_discrete(labels=feats_per_run_samples$sample) +
      theme_bw() + 
      theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
            axis.text.x=element_text(size=10), axis.text.y=element_text(angle=90),
            legend.text=element_text(size=10), legend.title=element_blank())
    # plotly doesnt support hjust, so calculate position to be in middle
    ggp = ggp + geom_text(aes(x=Run, y=nr_feats / 2, label=nr_feats), size=10, colour="white")
    p = ggplotly(ggp, width=600, height=vert_height) %>%
            layout(legend = list(orientation = 'h', x = 0, y = 1.1, xanchor='left', yanchor='bottom'))
    # Work around since plotly does not honor above legend.title=element_blank call
    p$x$layout$legend$title$text = ''
    htmlwidgets::saveWidget(p, glue('{feattype}plothtml/nrfeats__{chunk}.html'), selfcontained=F)
  
    # Feature count plot
  
    # precursor count per feat w quant plot
    nrprec = aggregate(precfeats_filt[[quantcol]], by=precfeats_filt[c(featcol, 'sample', filenamecol)], length)
    colnames(nrprec) = c(featcol, 'sample', filenamecol, 'Nr.precursors')
  
    # Box plots:
    plots = list(
      quant=list(feats=feats_filtered, col=quantcol),
      nrp=list(feats=nrprec, col='Nr.precursors')
    )
    for (plot in names(plots)) {
      data_to_plot = plots[[plot]]$feats
      col_to_plot = plots[[plot]]$col
      summary_stats <- data_to_plot[c('Run', col_to_plot)] %>%
        group_by(Run) %>%
        boxplot_stats(.data[[col_to_plot]])
      summary_stats$Run = as.character(summary_stats$Run)
      sumstats_samples = merge(summary_stats, inputfn[,c('Run', samplecol)], by='Run')
  
      # Plot using geom_crossbar, geom_errorbar
      ggp = ggplot(summary_stats, aes(x=Run, y=middle)) +
        geom_errorbar(aes(ymin = whisk_min, ymax = whisk_max), width=0) +
        geom_crossbar(aes(ymin = lower, ymax=upper), fill='white', linewidth=0.15) +
        coord_flip() + ylab(ptypes[[ptype]][2]) + theme_bw() + 
        ylab(col_to_plot) +
        scale_x_discrete(labels=sumstats_samples$sample) +
        theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
  	    axis.text=element_text(size=10), axis.text.y=element_text(angle=90)) +
        geom_text(aes(x=Run, y=(get('whisk_min') + get('whisk_max')) / 2,
          label=Run), position=position_nudge(x=0.1), size=3, colour="black")
      p = ggplotly(ggp, width=600, height=vert_height)
      htmlwidgets::saveWidget(p, glue('{feattype}plothtml/{plot}__{chunk}.html'), selfcontained=F)
    }

  for (nn_precursortable in nn_precursors) {
    precs_labeled = read.table(nn_precursortable, header=T, sep="\t", comment.char = "", quote = "")
    featcol = featmap[[feattype]][1]
    featfiltcol = featmap[[feattype]][2]
    quantcol = featmap[[feattype]][3]
    prec_featcols = precs_labeled[, c('Run', samplecol, featmap[[feattype]]) ]
    precfeats_filt = prec_featcols[prec_featcols[featfiltcol] < opt$conflvl, ]
    # Keep only one precursor line per feature per run
    feats = precfeats_filt %>%
      distinct(Run, .data[[featcol]], .keep_all=T)
  
    # Remove NA/0 for both precursor and feats table:
    precfeats_filt = precfeats_filt[precfeats_filt[[quantcol]] >= 0,]
    feats_filtered = feats[feats[[quantcol]] >= 0,]

    summary_stats <- feats_filtered[c('Run', quantcol)] %>%
      group_by(Run) %>%
      boxplot_stats(.data[[quantcol]])
    summary_stats$Run = as.character(summary_stats$Run)
    sumstats_samples = merge(summary_stats, inputfn[,c('Run', samplecol)], by='Run')

    # Plot using geom_crossbar, geom_errorbar
    ggp = ggplot(summary_stats, aes(x=Run, y=middle)) +
      geom_errorbar(aes(ymin = whisk_min, ymax = whisk_max), width=0) +
      geom_crossbar(aes(ymin = lower, ymax=upper), fill='white', linewidth=0.15) +
      coord_flip() + ylab(ptypes[[ptype]][2]) + theme_bw() + 
      ylab(quantcol) +
      scale_x_discrete(labels=sumstats_samples$sample) +
      theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(),
	    axis.text=element_text(size=10), axis.text.y=element_text(angle=90)) +
      geom_text(aes(x=Run, y=(get('whisk_min') + get('whisk_max')) / 2,
        label=Run), position=position_nudge(x=0.1), size=3, colour="black")
    p = ggplotly(ggp, width=600, height=vert_height)
    htmlwidgets::saveWidget(p, glue('{feattype}plothtml/nnquant__{chunk}.html'), selfcontained=F)
  }
  chunk = chunk + 1
  }

  # Text about overlap above protein/gene plots
  if (ncol(accumulated_feats) > 2) {
    accumulated_feats$rowsums = rowSums(accumulated_feats[,-1], na.rm=T)
    accumulated_feats_q$rowsums = rowSums(accumulated_feats_q[,-1], na.rm=T)
  } else {
    colnames(accumulated_feats) = c('feat', 'rowsums')
    colnames(accumulated_feats_q) = c('feat', 'rowsums')
  }
  feat_set_overlap = aggregate(feat~rowsums, accumulated_feats, length)
  feat_set_overlap_q = aggregate(feat~rowsums, accumulated_feats_q, length)
  write.table(feat_set_overlap, glue('{feattype}__overlap'), row.names=F, quote=F, sep='\t')

  if (nrow(inputfn) == max(feat_set_overlap$rowsums)) {
    all_run_overlap = feat_set_overlap$feat[feat_set_overlap$rowsums == nrow(inputfn)]
  } else {
    all_run_overlap = 0
  }
  if (nrow(inputfn) == max(feat_set_overlap_q$rowsums)) {
    overlap_q = feat_set_overlap_q$feat[feat_set_overlap_q$rowsums == nrow(inputfn)]
  } else {
    overlap_q = 0
  }
  writeLines(c(glue('Overlap for all runs: {all_run_overlap}'),
    glue('Overlap for all runs with quant: {overlap_q}'),
    glue('Total identified: {nrow(accumulated_feats)}')),
    glue('{feattype}plothtml/nrfeats__text.html'))
}

