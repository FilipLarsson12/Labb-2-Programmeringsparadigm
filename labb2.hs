-- Filip Larsson
module F2 where   

import Data.List
import GHC.Exts.Heap (GenClosure(prof))

-- 2.1
data MolType = Protein | DNA deriving (Eq, Show)
data MolSeq = MolSeq {name :: String, sequence :: String, molType :: MolType } deriving (Show, Eq)

listfordna = [0,1,2,3]
listforprotien = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

a = Profile{matrix=[[('A', 4), ('B', 0), ('C', 0), ('D', 0)], [('A', 1), ('B', 0), ('C', 1), ('D', 2)]], typ=DNA, sequenceAmount=4, nameofProfile="profile2"}
b = Profile{matrix=[[('A', 0), ('B', 1), ('C', 2), ('D', 1)], [('A', 2), ('B',2), ('C', 0), ('D', 0)]], typ=DNA, sequenceAmount=4, nameofProfile="profile1"}

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


nucleotides :: String
nucleotides = "ACGT"
aminoacids = sort "ARNDCEQGHILKMFPSTWYV"
makeProfileMatrix :: [MolSeq] -> [[(Char, Int)]]
makeProfileMatrix [] = error "Empty sequence list"
makeProfileMatrix sl = res
    where
        t = seqType (head sl)
        defaults =
            if (t == DNA) then
                zip nucleotides (replicate (length nucleotides) 0) -- Rad (i)
            else
                zip aminoacids (replicate (length aminoacids) 0) -- Rad (ii)
        strs = map seqSequence sl -- Rad (iii)
        tmp1 = map (map (\x -> ((head x), (length x))) . group . sort)
            (transpose strs) -- Rad (iv)
        equalFst a b = (fst a) == (fst b)
        res = map sort (map (\l -> unionBy equalFst l defaults) tmp1)


profileName :: Profile -> String
profileName (Profile _ _ _ nameofProfile) = nameofProfile

profileFrequency :: Profile -> Int -> Char -> Double
profileFrequency (Profile matrix _ sequenceAmount _) position letter = fromIntegral amount / fromIntegral sequenceAmount
    where 
        amount = matrixAmount (matrix !! position) letter

matrixAmount :: [(Char, Int)] -> Char -> Int 
matrixAmount [] _ = 0 
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
profileDistance profile1 profile2 = sumColumnDiffsOverRows profile1 profile2 (matrix profile1)

indexDiff :: Profile -> Profile -> Char -> Int -> Double
indexDiff profile1 profile2 char int =
    abs (profileFrequency profile1 int char - profileFrequency profile2 int char)

sumDiffsinColumn :: Profile -> Profile -> [(Char, Int)] -> Double
sumDiffsinColumn profile1 profile2 [] = 0
sumDiffsinColumn profile1 profile2 (h:t) =
    indexDiff profile1 profile2 (fst h) (snd h) + sumDiffsinColumn profile1 profile2 t

sumColumnDiffsOverRows :: Profile -> Profile -> [[(Char, Int)]] -> Double
sumColumnDiffsOverRows profile1 profile2 [] = 0
sumColumnDiffsOverRows profile1 profile2 (h:t) =
    sumDiffsinColumn profile1 profile2 h + sumColumnDiffsOverRows profile1 profile2 t

class Evol object where
    evolname :: object -> String 
    distance :: object -> object -> Double

    distancematrix :: [object] -> [(String, String, Double)]
    distancematrix [] = []
    distancematrix object = comparehead object 0 ++ distancematrix (tail object)


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


