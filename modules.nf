process DiaQuantificationReport {

  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple path(raws, arity: '1..*'), path('quants/*'), path(lib), path(fasta), val(diannparams), val(quantparams), val(enzyme), path(inputfn), val(normalize)
  
  output:
  tuple path('report.parquet'), path('*.tsv'), emit: report
  path('precursors_*.txt'), emit: precursors_split
  path('quantify_report.log'), emit: log

  script:
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)"}.join(' ') } \
    --lib $lib \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --use-quant \
    --quant-ori-names \
    --temp quants \
    ${quantparams.trim()} \
    ${diannparams.ms1acc ? "--mass-acc-ms1 ${diannparams.ms1acc}" : ''} \
    ${diannparams.ms2acc ? "--mass-acc ${diannparams.ms2acc}" : ''} \
    ${diannparams.window ? "--window $diannparams.window" : ''} \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --matrices \
    ${diannparams.ntermmetex ? '--met-excision' : ''} \
    ${diannparams.ntermac ? '--var-mod UniMod:1,42.010565,*n' : ''} \
    ${diannparams.idstonames ? '--ids-to-names' : ''} \
    --pg-level ${diannparams.pglvl} \
    --min-pr-charge ${diannparams.mincharge} \
    --max-pr-charge ${diannparams.maxcharge} \
    --min-pep-len ${diannparams.minpeplen} \
    --max-pep-len ${diannparams.maxpeplen} \
    --min-pr-mz ${diannparams.minmz} \
    --max-pr-mz ${diannparams.maxmz} \
    --min-fr-mz ${diannparams.minfrmz} \
    --max-fr-mz ${diannparams.maxfrmz} \
    ${diannparams.indiwin ? "--individual-windows" : ''} \
    ${diannparams.indiacc ? "--individual-mass-acc" : ''} \
    --qvalue ${diannparams.precfdr} \
    --matrix-qvalue ${diannparams.protfdr} \
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''} \
    ${normalize ? '': '--no-norm'} \
      | tee stdout.bak
    grep ERROR stdout.bak && exit 1
    parquet_to_tsv.py report.parquet $enzyme $inputfn ${normalize ? 'n' : 'nn'}

    mv report.log.txt quantify_report.log
  """
}
def read_header(info_fn) {
  def header = []
  def info = file(info_fn).eachLine { line, ix ->
      if (ix == 1) {
        header = line.tokenize('\t')
      }
  }
  return header
}


def create_info_map(info_fn, possible_params) {
  /* From possible params this parses a tab separated input file with
  a header (which has some of those params names.
  It returns a map with params and their values, which are defaults for
  those not set in the input file.
  */
  def info = file(info_fn).readLines().collect { it.tokenize('\t') }
  def header = info.pop()

  def params_not_header = possible_params - header

  def info_map = [:]
  def fpath
  def key_fn
  def tmp_fn
  info.findAll{ it[0][0] != '#' }.eachWithIndex { it, ix ->
    tmp_fn = [:]
    header.eachWithIndex{ hfield, hix ->
    tmp_fn[hfield] = it[hix]
    }
    fpath = file(tmp_fn.file_path)
    key_fn = ix //tmp_fn.file_path
    info_map[key_fn] = tmp_fn
    info_map[key_fn].id = ix
    info_map[key_fn].file_path = fpath
    info_map[key_fn].filename = "${info_map[key_fn].file_path.baseName}.${info_map[key_fn].file_path.extension}"
    params_not_header.each {
      info_map[key_fn][it] = params[it]
    }
  }
  return info_map
}


def identify_info_map(info_fn) {
  def expected_fields = ["file_path", "create_lib", "train_quantums", "quantfile"]
  def samples = create_info_map(info_fn, expected_fields)
  // Set all files to go to lib create if none specified (and no params.library passed)
  def new_samples = samples
  if (!samples.findAll { k,v -> v.create_lib as Integer }) {
    samples.each { k,v ->
      new_samples[k].create_lib = 1
    }
    samples = new_samples
  }
  // Set all files to go to train quantUMS if none specified
  if (!samples.findAll { k,v -> v.train_quantums as Integer }) {
    samples.each { k,v ->
      new_samples[k].train_quantums = 1
    }
    samples = new_samples
  }
  return samples
}


def listify(it) {
  /* This function is useful when needing a list even when having a single item
  - Single items in channels get unpacked from a list
  - Processes expect lists. Even though it would be fine
  without a list, for single-item-lists any special characters are not escaped by NF
  in the script, which leads to errors. See:
  https://github.com/nextflow-io/nextflow/discussions/4240
  */
  return it instanceof java.util.List ? it : [it]
}
