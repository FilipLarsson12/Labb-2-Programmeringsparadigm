module F2 where 

import Data.List
--import GHC.Exts.Heap (GenClosure(prof))

-- 2.1
data MolType = Protein | DNA deriving (Eq, Show)
data MolSeq = MolSeq {name :: String, sequence :: String, molType :: MolType } deriving (Show, Eq)


-- 2.2 
string2seq :: String -> String -> MolSeq
string2seq name seq
    | all (`elem` "ACGT") seq = MolSeq name seq DNA
    | otherwise = MolSeq name seq Protein

-- 2.3
seqName :: MolSeq -> String
seqName (MolSeq name _ _) = name

seqSequence :: MolSeq -> String
seqSequence (MolSeq _ seq _) = seq

seqLength :: MolSeq -> Int
seqLength = length . seqSequence

seqType :: MolSeq -> MolType
seqType (MolSeq _ _ moltype) = moltype

-- 2.4
seqDistance :: MolSeq -> MolSeq -> Double
seqDistance molseq1 molseq2
  | molType molseq1 /= molType molseq2 = error "Can't compare sequences of different molecular types"
  | molType molseq1 == DNA = dnaDistance (seqSequence molseq1) (seqSequence molseq2)
  | molType molseq1 == Protein = proteinDistance (seqSequence molseq1) (seqSequence molseq2)
  where
    dnaDistance s1 s2
      | alpha <= 0.74 = (-3/4) * log (1 - (4 * alpha / 3))
      | otherwise = 3.3
      where alpha = hammingDistance s1 s2 / fromIntegral (length s1)

    proteinDistance s1 s2
      | alpha <= 0.94 = (-19/20) * log (1 - (20 * alpha / 19))
      | otherwise = 3.7
      where alpha = hammingDistance s1 s2 / fromIntegral (length s1)

    hammingDistance s1 s2 = fromIntegral $ length $ filter (uncurry (/=)) $ zip s1 s2

-- 3.1 
data Profile = Profile {matrix :: [[(Char, Int)]], typ :: MolType, sequenceAmount :: Int, nameofProfile :: String} deriving (Show, Eq)

nucleotides = "ACGT"
aminoacids = sort "ARNDCEQGHILKMFPSTWYV"
makeProfileMatrix :: [MolSeq] -> [[(Char, Int)]]
makeProfileMatrix [] = error "Empty sequence list"
makeProfileMatrix sl = res
    where
        t = seqType (head sl)
        defaults =
            if (t == DNA) then
                zip nucleotides (replicate (length nucleotides) 0) ---- Rad (i) skapar lista med par. Varje par i listan har en bokstav och en nolla. Tex [('A', 0), ('C', 0), ('G', 0), ('T', 0)]
            else
                zip aminoacids (replicate (length aminoacids) 0) -- Rad (ii) skapar lista med par. Varje par i listan har en bokstav och en nolla.
        strs = map seqSequence sl -- Rad (iii) // skapar en lista som innehåller alla sekvenser från inlistan med Moltypes, tex: ['ACCT', 'CATT', 'TGCA'].
        tmp1 = map (map (\x -> ((head x), (length x))) . group . sort)
            (transpose strs) -- Rad (iv)
        -- Indata: ["ACCT", "CATT", "TGCA"]. 
        -- 1. tranpose ["ACCT", "CATT", "TGCA"] = ["ACT", "CAG", "CTC", "TTA"].
        -- 2. sort ["ACT", "CAG", "CTC", "TTA"] = ["ACT", "ACG", "CCT", "ATT"].
        -- 3. group ["ACT", "ACG", "CCT", "ATT"] = [["A", "C", "T"], ["A", "C", "G"], ["CC", "T"], ["A", "TT"]]
        -- 4. map (\x -> ((head x), (length x))) [[('A', 1), ('C', 1), ('T', 1)], [('A', 1), ('C', 1), ('G', 1)], [('C', 2), ('T', 1)], [('A', 1), ('T', 2)]]
        -- Hela instruktionen: Först transponerar vi listan så vi har en lista med första bokstäverna i sekvensen, andra bokstäverna,
        -- osv osv. Sedan Sorterar vi första bokstäverna efter bokstavordning, sorterar andra bokstäverna efter bokstavsordning, 
        -- tredje osv osv. Sedan grupperar vi första bokstäverna, andra bokstäverna, tredje osv osv.
        -- Tex blir "AAACCT" = "AAA", "CC", "T".
        -- Sedan så för varje sådan lista med bokstäver så skapar vi en tuple med head av listan samt listans längd. På så sätt
        -- får vi hela matrixen. "AAA", "CC", "T" = ('A', 3), ('C', 2), ('T', 1)
        equalFst a b = (fst a) == (fst b)
        res = map sort (map (\l -> unionBy equalFst l defaults) tmp1)


profileName :: Profile -> String
profileName (Profile _ _ _ nameofProfile) = nameofProfile

profileFrequency :: Profile -> Int -> Char -> Double
profileFrequency (Profile matrix _ sequenceAmount _) position letter = fromIntegral amount / fromIntegral sequenceAmount
    where 
        amount = matrixAmount (matrix !! position) letter 

matrixAmount :: [(Char, Int)] -> Char -> Int 
matrixAmount (h:t) (letter)
    | fst h == letter = snd h 
    | otherwise = matrixAmount t letter 

molseqs2profile :: String -> [MolSeq] -> Profile
molseqs2profile s molseqlist = Profile matrix typ sequenceAmount nameofProfile
    where
        matrix = makeProfileMatrix molseqlist
        typ = seqType (head molseqlist)
        sequenceAmount = length molseqlist
        nameofProfile = s

profileDistance :: Profile -> Profile -> Double
profileDistance profile1 profile2 = sumMatrices profile1 profile2 matrix1 matrix2
    where 
        matrix1 = matrix profile1
        matrix2 = matrix profile2

sumMatrices :: Profile -> Profile -> [[(Char, Int)]] -> [[(Char, Int)]] -> Double 
sumMatrices _ _ [] []  = 0
sumMatrices profile1 profile2 matrix1 matrix2 = columnDiffResult + sumMatrices profile1 profile2 (tail matrix1) (tail matrix2)
    where 
        columnDiffResult = columnDiff (head matrix1) (head matrix2) profile1 profile2

columnDiff :: [(Char, Int)] -> [(Char, Int)] -> Profile -> Profile -> Double 
columnDiff [] [] _ _ = 0
columnDiff column1 column2 profile1 profile2 = abs(m1 - m2)  + columnDiff (tail column1) (tail column2) profile1 profile2 
    where 
        m1 = fromIntegral (snd (head column1)) / fromIntegral (sequenceAmount profile1)
        m2 = fromIntegral (snd (head column2)) / fromIntegral (sequenceAmount profile2)




class Evol object where
    evolname :: object -> String 
    distance :: object -> object -> Double

    distanceMatrix :: [object] -> [(String, String, Double)]
    distanceMatrix [] = []
    distanceMatrix object = comparehead object 0 ++ distanceMatrix (tail object)


    comparehead :: [object] -> Int -> [(String, String, Double)]
    comparehead object number 
        | number < length object = (evolname headname, evolname tailname , distance headname tailname) : comparehead object (number+1)
        | otherwise = [] 
        where 
            headname = head object 
            tailname = object !! number 


instance Evol MolSeq where
    evolname = seqName
    distance = seqDistance

instance Evol Profile where 
    evolname = nameofProfile
    distance = profileDistance


