// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;


SEXP generic_logical_subset( SEXP xin , LogicalVector w);


template <int RTYPE>
Vector<RTYPE> generic_logical_subset_impl( Vector<RTYPE> xin, LogicalVector w){
  return xin[w] ;
}

SEXP generic_logical_subset( SEXP xin , LogicalVector w){
    RCPP_RETURN_VECTOR(generic_logical_subset_impl, xin, w) ;
}

// [[Rcpp::export]]
IntegerVector naomitwhy(DataFrame df, Function is_na_generic) {
  int m = df.nrow();
  int n = df.ncol();
  CharacterVector df_names = df.names();


  LogicalVector omit = LogicalVector(m);

  List why_omit(n);
  why_omit.names() = df_names;
  LogicalVector why_omit_idx(n);

  bool anyomit = false;

  for (int j =0; j < n; j++) {
    std::string nm = as<std::string>(df_names(j));

    LogicalVector isna = is_na_generic(df(j));

    SEXP d  = isna.attr("dim");

    if(!Rf_isNull(d) && Rf_length(d) == 2  && INTEGER(d)[1] > 1){
      int dcols = INTEGER(d)[1];


      LogicalVector na2(m);

      for(int k = 0, ii=0; k < dcols; k++)
      for(int i = 0; i < m; i++, ii++){
        if(isna[ii]) na2[i] = true;
      }

      isna = na2;
    }

    int ii = 0;
    IntegerVector why_omit_j(m);

    for (int i = 0; i < m; i++){

      if(isna[i]){
        if(!omit[i]){
          why_omit_j[ii++] = i + 1;
        }

        omit[i] = true;
      }
    }

    if(ii > 0){
      why_omit[j] = why_omit_j[seq(0,ii-1)];
      why_omit_idx[j] = true;
      anyomit = true;
    }
  }
  if(!anyomit){ return(NULL); }

  IntegerVector omit_idx = seq_len(m);
  omit_idx = generic_logical_subset(omit_idx, omit);
  omit_idx.attr("names") = generic_logical_subset(df.attr("row.names"), omit);
  omit_idx.attr("why_omit") = why_omit[why_omit_idx];
  omit_idx.attr("class") = CharacterVector::create("omit", "detailed");



  return(omit_idx);
}

// # NJF 10/18
// # Silly microbenchmark to make sure I didn't make it slower
// # df <- expand.grid(x=c(1:100, NA), y=c(1:5, NA), z=c(1:8, NA), q=c(NA,2:5))
//
// # microbenchmark(stock=na.omit(df), hack1=na.omit_detailed.data.frame(df))
// ## Unit: milliseconds
// ## expr      min       lq     mean   median       uq        max neval
// ## stock 6.114132 6.184318 7.744881 6.232744 6.961491 101.823530   100
// ## hack1 5.360638 5.480531 6.525075 5.694078 7.752104   9.323943   100

// > microbenchmark(
//     +   out1 = estimatr:::naomitwhy(df),
//     +   out2 = estimatr:::na.omit_detailed.data.frame(df)
//     + )
//     Unit: milliseconds
//   expr      min        lq       mean    median       uq       max neval
//   out1 2.830644 2.9078700 3.88764381 2.9900055 3.655617 67.368491   100
// out2 3.463163 3.7196865 5.70407818 4.5768560 4.748749 70.595884   100
// > df <- na.omit(df)
//   > nrow(df)
//   [1] 16000
// > microbenchmark(
//     +   out1 = estimatr:::naomitwhy(df),
//     +   out2 = estimatr:::na.omit_detailed.data.frame(df)
//     + )
//     Unit: microseconds
//   expr      min        lq       mean    median        uq       max neval
//   out1  152.816  159.3180  214.70463  171.3470  185.2210  1023.877   100
// out2 1200.924 1227.9935 2193.33799 1271.2825 2040.2235 64933.756   100
// >
