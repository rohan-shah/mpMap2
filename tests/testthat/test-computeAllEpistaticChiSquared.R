context("computeAllEpistaticChiSquared")
test_that("Check that computeAllEpistaticChiSquared agrees with the mpMap version, for a 4-way design",
{
	map <- qtl::sim.map(len = 5, n.mar = 6, anchor.tel = TRUE, include.x=FALSE, eq.spacing=TRUE)
	pedigree <- mpMap::sim.mpped(nfounders = 4, nfunnels = 3, nperfam = 500, nssdgen = 6, iripgen = 0, seed = 1)
	mpMapObject <- mpMap::sim.mpcross(map = map, pedigree = pedigree)
	mpMapObject$map <- map
	capture.output(mpMapObject <- mpMap::mpprob(mpMapObject, step = 0, mrkpos = TRUE, mapfx = "haldane", est = FALSE, ibd = FALSE))
	unlink("tmp.ril.csv")
	unlink("tmp.founder.csv")

	capture.output(csfx03 <- inline::cxxfunction(signature(x="numeric", y="numeric", nf="numeric"), plugin="RcppArmadillo", body="
		arma::mat xv = Rcpp::as<arma::mat>(x);
		arma::mat yv = Rcpp::as<arma::mat>(y);
		int nfou=Rcpp::as<int>(nf);
		int nm1=yv.n_cols/nfou, nry=yv.n_rows, nm2=xv.n_cols/nfou, i, j, k, m;
		arma::mat z = arma::randu<arma::mat>(nm1, nm2);
		arma::vec t = arma::randu<arma::vec>(nm1);
		arma::vec exp = arma::randu<arma::vec>(nfou*nfou);
		arma::vec obs = arma::randu<arma::vec>(nfou*nfou);

		for (j=0; j<nm2; j++) 
		{
			t.zeros();
			for (i=0; i<=j; i++)
			{
				for (k=0; k<nfou; k++)
				for (m=0; m<nfou; m++) 
				{
					obs(k*nfou+m) = sum(yv.col(j*nfou+k)%xv.col(i*nfou+m));
					exp(k*nfou+m) = sum(yv.col(j*nfou+k))*sum(xv.col(i*nfou+m))/nry;
				}
				t(i) = sum((abs(obs-exp)-.5)%(abs(obs-exp)-.5)/exp);
			}
			z.col(j) = t;
		}
		return Rcpp::List::create(Rcpp::Named(\"ts\")=z);
		"
	))
	mpMapProbabilities <- do.call(cbind, mpMapObject$prob)
	mpMapChiSquared <- csfx03(mpMapProbabilities, mpMapProbabilities, 4)$ts
	gdata::lowerTriangle(mpMapChiSquared) <- gdata::lowerTriangle(t(mpMapChiSquared))

	suppressWarnings(mpMap2Object <- mpMap2::fromMpMap(mpMapObject, selfing = "infinite", fixCodingErrors = FALSE))
	mpMap2Object <- new("mpcrossMapped", mpMap2Object, map = map)
	mpMap2Probabilities <- computeGenotypeProbabilities(mpMap2Object)
	mpMap2ChiSquared <- .Call("computeAllEpistaticChiSquared", mpMap2Probabilities@geneticData[[1]]@probabilities, 4, TRUE, PACKAGE="mpMap2")

	expect_true(all(abs(mpMap2ChiSquared - mpMapChiSquared) / mpMapChiSquared < 0.01))
	expect_true(all(abs(mpMap2ChiSquared - mpMapChiSquared) / mpMap2ChiSquared < 0.01))
})
