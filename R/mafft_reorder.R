#' MAFFT Alignment with Reorder
#'
#' This function adds query sequences to a reference alignment and realigns all the sequences, allowing
#' reordering of the sequences
#'
#' @param sequences The query sequences to be added to the reference alignment in fasta format
#' @param ref_align The reference alignment to add the query sequences to in fasta format
#' @return An alignment containing the query sequences and reference sequences
#' @export
mafft_reorder<-function(sequences, ref_align){
  OS<-.Platform$OS.type

  if (OS == "unix")
    exec <- "/usr/local/bin/mafft"
  if (OS == "windows")
    exec <- "C:\\Windows\\System32\\wsl.exe mafft"

  sequences <- as.list(sequences)
  fns <- vector(length = 4)
  for (i in seq_along(fns)) fns[i] <- tempfile(pattern = "mafft",
                                               tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])


  ref_align <- as.list(ref_align)


  ips::write.fas(sequences, fns[1])
  ips::write.fas(ref_align, fns[2])
  call.mafft <- paste(exec, " --add ", fns[1], " --reorder ", fns[2], " > ",
                      fns[3], sep = "")

  if (OS == "unix") {
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- (file.info(fns[3])$size > 1)
    if (res != 0) {
      res <- ape::read.dna(fns[3], format = "fasta")
    }
  }else{
    res <- system(call.mafft, intern = TRUE, ignore.stderr = FALSE)
    if (length(grep("error|ERROR", res))) {
      res <- 0
    } else {
      res <- ape::read.dna(fns[3], format = "fasta")
    }
  }
  unlink(fns[file.exists(fns)])

  res<-as.character(res)

  return(res)

}









