context("Test reorderPedigree function")
test_that("Test reorderPedigree function using permutations",
	{
		for(i in 1:10)
		{
			pedigree <- eightParentPedigreeRandomFunnels(initialPopulationSize = 100, selfingGenerations = 1, nSeeds = 1, intercrossingGenerations = 0)
			pedigree <- as(pedigree, "pedigree")
			permutation <- sample(1:length(pedigree@lineNames))
			newLineNames <- pedigree@lineNames[permutation]
			tmp <- pedigree@mother[match(newLineNames, pedigree@lineNames)]
			tmp[tmp == 0] <- NA
			newMotherLineNames <- pedigree@lineNames[tmp]
			newMother <- match(newMotherLineNames, newLineNames)
			newMother[is.na(newMother)] <- 0L
		
			tmp <- pedigree@father[match(newLineNames, pedigree@lineNames)]
			tmp[tmp == 0] <- NA
			newFatherLineNames <- pedigree@lineNames[tmp]

			newFather <- match(newFatherLineNames, newLineNames)
			newFather[is.na(newFather)] <- 0L
			expect_error(new("pedigree", lineNames = newLineNames, mother = newMother, father = newFather, selfing = pedigree@selfing, warnImproperFunnels = pedigree@warnImproperFunnels))

			reordered <- mpMap2:::reorderPedigree(mother = newMother, father = newFather, lineNames = newLineNames, selfing = pedigree@selfing, warnImproperFunnels = pedigree@warnImproperFunnels)
			if(!is.null(reordered))
			{
				#Now check that the mother and father of every line are the same
				reorderedMotherLineNames <- reordered@lineNames[reordered@mother]
				reorderedFatherLineNames <- reordered@lineNames[reordered@father]
				expect_identical(reorderedMotherLineNames[match(pedigree@lineNames[-(1:8)], reordered@lineNames)-8], pedigree@lineNames[pedigree@mother])
				expect_identical(reorderedFatherLineNames[match(pedigree@lineNames[-(1:8)], reordered@lineNames)-8], pedigree@lineNames[pedigree@father])
			}
		}
	})
