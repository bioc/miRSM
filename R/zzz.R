.onAttach <- function(lib, pkg) { 
    
    version <- packageDescription("miRSM", fields = "Version")
    packageStartupMessage( "Citation: [1] Zhang J, Liu L, Xu T, Zhang W, Zhao C, Li S, Li J, Rao N, Le TD,","\n",
      "miRSM: an R package to infer and analyze miRNA sponge modules in heterogeneous data,","\n",
      "RNA Biol. 2021, 18(12):2308-2320.","\n",
      "[2] Zhang J, Wei X, Zhao C, Yang H,","\n",
      "Protocol to infer and analyze miRNA sponge modules in heterogeneous data using miRSM 2.0,","\n",
      "STAR Protoc. 2024, 5(4):103317.","\n",
      "BibTex: enter 'toBibtex(citation(\"miRSM\"))","\n\n",
      "Homepage: https://github.com/zhangjunpeng411/miRSM","\n\n",
      "miRSM Package Version ", version, "\n")
}
