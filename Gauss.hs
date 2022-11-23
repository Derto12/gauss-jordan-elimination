{-# LANGUAGE MultiWayIf #-}

module Gauss where

import Data.List (intercalate)
import Data.Maybe

-- GElem and GElem related functions
data GElem = G Rational Bool Bool
 deriving (Eq)

instance Ord GElem where
    (>) (G a1 _ _) (G a2 _ _) = a1 > a2
    (<) (G a1 _ _) (G a2 _ _) = a1 < a2
    (>=) (G a1 _ _) (G a2 _ _) = a1 >= a2
    (<=) (G a1 _ _) (G a2 _ _) = a1 <= a2

instance Num GElem where
    (+) (G a1 marked mNeigh) (G a2 _ _) = (G (a1 + a2) marked mNeigh)
    (-) (G a1 marked mNeigh) (G a2 _ _) = (G (a1 - a2) marked mNeigh)
    (*) (G a1 marked mNeigh) (G a2 _ _) = (G (a1 * a2) marked mNeigh)
    abs (G a marked mNeigh) = (G (abs a) marked mNeigh)
    signum (G a marked mNeigh) = (G (signum a) marked mNeigh)
    fromInteger i = (G (fromInteger i) False False)

instance Fractional GElem where
    (/) (G a marked mNeigh) (G b _ _) = (G (a/b) marked mNeigh)
    fromRational r = (G (fromRational r) False False)

instance Show GElem where
    show (G a marked mNeigh)
     | marked = "<("++strA++")>"
     | otherwise = "("++strA++")"
      where
      strA = 
        if drop c str == " % 1" 
        then take c str
        else show a
         where
          str = show a
          c = length str - 4
    showList (xs) = ((intercalate " " $ map show xs) ++)

toGelem :: Rational -> GElem
toGelem a = G a False False

toGList :: [Rational] -> [GElem]
toGList as = map toGelem as

toMRowList :: [Rational] -> [GElem]
toMRowList as = map (\a -> (G a False True)) as

toGNestedList :: [[Rational]] -> [[GElem]]
toGNestedList as = map toGList as

toRational :: GElem -> Rational
toRational (G a False False) = a

printTable :: Int -> [[GElem]] -> IO () -- where the Int represents the number of columns on the LHS
printTable i as = putStrLn $ unlines $ map (\a -> (showRow $ take i a) ++ "|    " ++ (showRow $ drop i a)) as
 where
 showRow as = flatten $ map (\a -> show a ++ replicate (8 - (length $ show a)) ' ') as

-- Other helper functions
flatten :: [[a]] -> [a]
flatten as = foldr (++) [] as

changeElem :: Int -> a -> [a] -> [a]
changeElem i a as = as1 ++ a:as2
  where
  (as1, _:as2) = splitAt i as

addElemAt :: Int -> a -> [a] -> [a]
addElemAt i a as = as1 ++ a:as2
  where
  (as1, as2) = splitAt i as

removeElemAt :: Int -> [a] -> [a]
removeElemAt i as = as1 ++ as2
  where
  (as1, _:as2) = splitAt i as

findWhere :: (a -> Bool) -> [a] -> Maybe (Int, a)
findWhere _ [] = error "empty list"
findWhere f ls = findHelper f 0 ls
    where
    findHelper ::(a -> Bool) ->  Int -> [a] -> Maybe (Int, a)
    findHelper _ _ [] = Nothing
    findHelper f currInd (a:as)
     | f a = Just (currInd, a)
     | otherwise = findHelper f (currInd+1) as

findWhereMulti :: (Ord a) => (a -> Bool) -> [[a]] -> Maybe ((Int, Int), a)
findWhereMulti f as
 | isJust maybeA = Just (newInd, val)
 | otherwise = Nothing
    where
    maybeA = findWhere f (flatten as)
    (ind, val) = fromJust maybeA
    newInd = (div ind (length (as !! 0)), mod ind (length (as !! 0)))

hasSameLengthLists :: [[a]] -> Bool
hasSameLengthLists (a:as) = hasHelper a as
    where
    hasHelper :: [a] -> [[a]] -> Bool
    hasHelper _ [] = True
    hasHelper a (l:ls)
     | length a == length l = hasHelper a ls
     | otherwise = False


-- Gauss-Jordan elimination related functions
selectNext :: [[GElem]] -> Maybe ((Int, Int), GElem)
selectNext as = findWhereMulti (\(G a _ mNeigh) -> not mNeigh && a /= 0) as

markSelected :: (Int, Int) -> GElem -> [[GElem]] -> [[GElem]]
markSelected (i, j) (G v _ _) as = changeElem i (changeElem j (G v True True) (as!!i)) as

markRow :: Int -> [[GElem]] -> [[GElem]]
markRow j as = map (\a -> markHelper j (a!!j) a) as
    where
    markHelper :: Int -> GElem -> [GElem] -> [GElem]
    markHelper i (G a marked _) as = changeElem i (G a marked True) as

markColumn :: Int -> [[GElem]] -> [[GElem]]
markColumn i as = changeElem i a as
    where
    a = map (\(G a marked _) -> (G a marked True)) (as!!i)

markNeighbours :: (Int, Int) -> GElem -> [[GElem]] -> [[GElem]]
markNeighbours (i, j) g as = markColumn i (markRow j (markSelected (i, j) g as))

divideSelectedRow :: Int -> GElem -> [[GElem]] -> [[GElem]]
divideSelectedRow i a as = changeElem i (map (/a) (as!!i)) as

decrOtherRows :: (Int,Int) -> GElem -> [[GElem]] -> [[GElem]]
decrOtherRows (i,j) g as = addElemAt i gs (decrHelper (removeElemAt i as)) 
    where
    gs = as!!i
    decrHelper :: [[GElem]] -> [[GElem]]
    decrHelper ls = map mapHelper ls
     where
     mapHelper l = zipWith zipHelper l gs
      where
      zipHelper e1 e2 = e1 - ((l!!j) / g) * e2

prepareTable :: Int -> [[GElem]] -> [[GElem]]
prepareTable i as = prepareHelper i (length (as!!0)) as
 where
 prepareHelper :: Int -> Int -> [[GElem]] -> [[GElem]]
 prepareHelper i n as = markRow (n - i) as

executeSteps :: [[GElem]] -> [[GElem]]
executeSteps as
 | isJust maybeFound = executeSteps (
    decrOtherRows (i,j) (G 1 True True) (
        divideSelectedRow i g (
            markNeighbours (i,j) g (
                markSelected (i,j) g as)
    )))
 | otherwise = as
    where
    maybeFound = selectNext as
    ((i, j), g) = fromJust maybeFound

-- gauss table with no right hand side
gaussNRH :: [[Rational]] -> [[GElem]]
gaussNRH [] = error "empty list"
gaussNRH as
 | hasSameLengthLists as = executeSteps $ toGNestedList as
 | otherwise = error "varying list lengths"

-- gauss table with vector at the right hand side
gaussVRH :: [[Rational]] -> [Rational] -> [[GElem]]
gaussVRH as bs
 | as == [] || bs == [] = error "empty list"
 | hasSameLengthLists as && length as == length bs = executeSteps (
    zipWith (\a b -> a++[b]) (toGNestedList as) (toMRowList bs))
 | otherwise = error "varying list lengths"

-- gauss table with matrix at the right hand side
gaussMRH :: [[Rational]] -> [[Rational]] -> [[GElem]]
gaussMRH as bs
 | as == [] || bs == [] = error "empty list"
 | hasSameLengthLists as 
    && hasSameLengthLists bs 
    && length as == length bs = executeSteps (
        zipWith (\a b -> a++b) (toGNestedList as) (map toMRowList bs))
 | otherwise = error "varying list lengths"

-- identity matrix at the right hand side
-- gaussIMRH:: [[Rational]] -> [[Rational]] -> [[GElem]]