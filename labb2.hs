module F2 where

-- Uppgift 1
data MolSeq = MolSeq String String String deriving Show


-- Uppgift 2    

checkString :: String -> String
checkString [] = "DNA"
checkString (x:xs) 
        | x `elem` "ACGT" = checkString xs
        | otherwise = "Protein"

string2seq :: String -> String -> MolSeq
string2seq str1 str2
    | checkString str2 == "DNA" = MolSeq str1 str2 "DNA"
    | otherwise = MolSeq str1 str2 "Protein"
    

