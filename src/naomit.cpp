// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppEigen.h>
using namespace Rcpp;


template <int RTYPE>
Vector<RTYPE> generic_logical_subset_impl( Vector<RTYPE> xin, LogicalVector w){
  return xin[w] ;
}

SEXP generic_logical_subset( SEXP xin , LogicalVector w){
  RCPP_RETURN_VECTOR(generic_logical_subset_impl, xin, w) ;
}

// [[Rcpp::export]]
DataFrame naomitwhy(DataFrame df, LogicalMatrix isna, Function recursive_subset) {
  int m = df.nrow();
  int n = df.ncol();
  int N = isna.ncol();


  CharacterVector df_names = df.names();

  IntegerVector na_to_col_map(n);
  if(N == n){
    std::fill(na_to_col_map.begin(), na_to_col_map.end(), 1);
  }
  else {
    Function dim("dim");

    for(int i = 0; i < n; i++){
      SEXP dfi = df[i];
      if(Rf_isVectorAtomic(dfi) && LENGTH(dfi) == m){
        na_to_col_map[i] = 1;
      } else {
        SEXP nc = dim(dfi);
        na_to_col_map[i] = Rf_isNull(nc) ? 1 : INTEGER(nc)[1];
      }
    }
  }

  LogicalVector omit = LogicalVector(m);

  int omit_count = 0, omit_f = m, omit_l = 0;

  List why_omit(n);
  why_omit.names() = df_names;
  LogicalVector why_omit_idx(n);


  for (int j = 0, ii = 0; j < n; j++) {

    std::vector<int> why_omit_j;

    for (int j_sub = na_to_col_map[j]; j_sub; j_sub--){
      for (int i = 0; i < m; i++, ii++){

        if(isna[ii]){
          if(!omit[i]){
            why_omit_j.push_back(i + 1);
          }

          omit[i] = true;
        }
      }
    }

    if(why_omit_j.size() > 0){
      if(na_to_col_map[j] > 1){
        std::sort(why_omit_j.begin(), why_omit_j.end());
      }
      why_omit[j] = wrap(why_omit_j);
      why_omit_idx[j] = true;
      omit_f = std::min(omit_f, why_omit_j.front());
      omit_l = std::max(omit_l, why_omit_j.back());
      omit_count += why_omit_j.size();
    }
  }
  if(omit_count == 0){ return(df); }

  // Rcout << "after\n" << omit_count << "\n";

  IntegerVector omit_idx = IntegerVector(omit_count);
  for(int i = omit_f-1, ii=0; i < omit_l; i++){
    if(omit[i]) omit_idx[ii++] = i+1;
  }

  CharacterVector rownames = df.attr("row.names");
  omit_idx.attr("names") = rownames[omit];

  omit_idx.attr("why_omit") = why_omit[why_omit_idx];
  omit_idx.attr("class") = CharacterVector::create("omit", "detailed");
  //omit_idx.attr("tokeep") = !omit;

  omit = !omit;

  List out(n);

  for(int i = 0; i < n; i++){
    SEXP dfi = df(i);
    if(Rf_isVectorAtomic(dfi) && LENGTH(dfi) == m){
      out[i] = generic_logical_subset(dfi, omit);
    } else {
      out[i] = recursive_subset(dfi, omit);
    }

  }

  out.names() = df_names;
  out.attr("row.names") = rownames[omit];
  out.attr("na.action") = omit_idx;
  out.attr("class") = df.attr("class");

  return(out);
}


  // require(microbenchmark)
  // df <- expand.grid(x=c(1:100, NA), y=c(1:5, NA), z=c(1:8, NA), q=c(NA,2:5))
  // df2 <- na.omit(df)
  // microbenchmark(stock=na.omit(df), ours=estimatr:::na.omit_detailed.data.frame(df))
  // microbenchmark(stock=na.omit(df2), ours=estimatr:::na.omit_detailed.data.frame(df2), unit="ms")
  //
  // df <- rbind(df, df2, df)
  // df2 <- rbind(df2, df2, df2)
  //
  // microbenchmark(stock=na.omit(df), ours=estimatr:::na.omit_detailed.data.frame(df), unit="ms")
  // microbenchmark(stock=na.omit(df2), ours=estimatr:::na.omit_detailed.data.frame(df2), unit="ms")

  // df <- cbind(df, df,df)
  // df2 <- cbind(df2, df2, df2)
  // microbenchmark(stock=na.omit(df), ours=estimatr:::na.omit_detailed.data.frame(df), unit="ms")
  // microbenchmark(stock=na.omit(df2), ours=estimatr:::na.omit_detailed.data.frame(df2), unit="ms")

  // sleep[c("sleep", "foo")] = list(sleep, matrix(1:40, 20))
  // sleep[cbind(c(1,5,9), c(2,1,3))] <- NA
  // sleep$sleep[cbind(1+c(1,5,9), c(2,1,3))] <- NA
  // sleep$foo[12,1] <- NA
  // attributes(estimatr:::na.omit_detailed.data.frame(sleep))
