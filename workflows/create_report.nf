process extractThermoScans {
  /* Used in DIA for Thermo, so we dont need to convert to mzML
  just to get the nr of scans */
  tag 'scanheadsman'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  path(raws)

  output:
  path('nrscans')

  script:
  """
  exitcode=0
  ${raws.collect { "timeout --preserve-status 5s dotnet /scanheadsman/ScanHeadsman.dll \$(realpath $it) > tmpfn || exitcode=\$? \
    && if [[ \$exitcode != 0 && \$exitcode != 143 ]] ; then exit \$exitcode ;fi \
    && echo ${it.baseName}\$'\t'\$(grep 'Processing scan [0-9]' tmpfn | sed 's/.* of //') >> nrscans"}.join('\n')}
  """
}


process getBrukerScanNumbers {
  tag 'sqlite'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  path(raws)

  output:
  path('nrscans')

  script:
  """
  ${raws.collect { "echo ${it.baseName}\$'\t'\$(sqlite3 ${it}/analysis.tdf 'SELECT COUNT(*) FROM Frames') >> nrscans"}.join(' && ')}
  """
}


process precursorPlot {
  tag 'ddamsproteomics'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple path('filescans????'), path(precursors), path(nonnorm_prec), path(inputfn), val(conflvl)
  
  output:
  tuple path('precursorplothtml'), path('*_qc.txt'), path('*__overlap'), path('genesplothtml'), path('proteinsplothtml')
  
  script:
  nonnorm_parsed = nonnorm_prec[0].name == 'NO__FILE' ? '' : "--non-norm-precursors $nonnorm_prec"
  // FIXME error if not finding these columns!
  """
  cat filescans* > concat_filescans
  mkdir -p precursorplothtml genesplothtml proteinsplothtml
  precursor_qc.R --precursors $precursors ${nonnorm_parsed} --inputfn $inputfn --conflvl $conflvl
  """
}


process summaryReport {
cache false

  tag 'ddamsproteomics'
  container Containers.containers[task.tag][workflow.containerEngine]


  input:
  tuple path('precplothtml'), path(summaries), path(feat_overlaps), path('genesplots'), path('proteinsplots')
  
  output:
  tuple path('report_groovy_template.html'), path('report_groovy_template_light.html'), path('libs.js')
  
  script:
  """
  # xargs removes trailing whitespace
  report_tables.py --version "${workflow.manifest.version}" --doi "${workflow.manifest.doi}" \
      --templatedir "$baseDir/assets"
  """
}


workflow QC_REPORT {
  take:
  raws_ftypes
  inputfn
  precursors_split
  non_norm_precursors_split
  proteinfdr
  
  main:
  raws_ftypes
  | filter { it[1] == 'thermo' }
  | map { it[0] }
  | collate(10) // split raw reading in batches of 10
  | filter { it.size() > 0 }
  | extractThermoScans
  
  raws_ftypes
  | filter { it[1] == 'bruker' }
  | map { it[0] }
  | collate(10) // split raw reading in batches of 10
  | filter { it.size() > 0 }
  | getBrukerScanNumbers
  | mix(extractThermoScans.out)
  | toList
  | toList
  | combine(precursors_split | toSortedList)  // FIXME need to order and line up w non norm
  | combine(non_norm_precursors_split | toSortedList)
  | combine(inputfn)
  | map { it + [proteinfdr] }
  | precursorPlot
  | summaryReport
  
  emit:
  summaryReport.out
}
