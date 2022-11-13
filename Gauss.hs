{-# LANGUAGE MultiWayIf #-}

module Gauss where

import Data.List (intercalate)

data GElem = G Rational Bool Bool
 deriving (Eq)

instance Ord GElem where
    (>) (G a1 _ _) (G a2 _ _) = a1 > a2
    (<) (G a1 _ _) (G a2 _ _) = a1 < a2
    (>=) (G a1 _ _) (G a2 _ _) = a1 >= a2
    (<=) (G a1 _ _) (G a2 _ _) = a1 <= a2

instance Num GElem where
    (+) (G a1 marked mRow) (G a2 _ _) = (G (a1 + a2) marked mRow)
    (-) (G a1 marked mRow) (G a2 _ _) = (G (a1 - a2) marked mRow)
    (*) (G a1 marked mRow) (G a2 _ _) = (G (a1 * a2) marked mRow)
    abs (G a marked mRow) = (G (abs a) marked mRow)
    signum (G a marked mRow) = (G (signum a) marked mRow)
    fromInteger i = (G (fromInteger i) False False)

instance Fractional GElem where
    (/) (G a marked mRow) (G b _ _) = (G (a/b) marked mRow)
    fromRational r = (G (fromRational r) False False)

instance Show GElem where
    show (G a marked mRow)
     | marked = "<<("++show a++")>>"
     | otherwise = "("++show a++")"
    showList (xs) = ((intercalate " " $ map show xs) ++)

-- TODO - better visualize data
printGauss :: [[GElem]] -> IO () 
printGauss as = putStrLn (unlines (map (\a -> show a) as))


toGelem :: Rational -> GElem
toGelem a = G a False False

toGList :: [Rational] -> [GElem]
toGList as = map toGelem as

toRational :: GElem -> Rational
toRational (G a False False) = a

joinLists :: [[a]] -> [a]
joinLists as = foldr (++) [] as

-- FIX - minWhere returns the minimum value, but the closest number to 1 would be more ideal
minWhere :: Ord a => (a -> Bool) -> [a] -> (Bool, Int, a)
minWhere _ [] = error "empty list"
minWhere f l@(a:as) = minHelper f False 0 a 0 l
    where
    minHelper :: Ord a => (a -> Bool) -> Bool -> Int -> a -> Int -> [a] -> (Bool, Int, a)
    minHelper _ found minInd minVal _ [] = (found, minInd, minVal)
    minHelper f found minInd minVal currInd (a:as)
     | found && f a = if
        | a < minVal -> minHelper f found currInd a (currInd+1) as
        | otherwise -> minHelper f found minInd minVal (currInd+1) as
     | not found && f a = minHelper f True currInd a (currInd+1) as
     | otherwise = minHelper f found minInd minVal (currInd+1) as

minWhereMulti :: (Ord a) => (a -> Bool) -> [[a]] -> (Bool, (Int, Int), a)
minWhereMulti f as = (found, newInd, val)
    where
    (found, ind, val) = minWhere f (joinLists as)
    newInd = (div ind (length (as !! 0)), mod ind (length (as !! 0)))

gaussSelector :: [[GElem]] -> (Bool, (Int, Int), GElem)
gaussSelector as = (found, (x, y), a)
 where (found, (x, y), a) = minWhereMulti (\(G a _ mRow) -> not mRow && a /= 0) as

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

uniteElems :: (a -> a -> a) -> [a] -> [a] -> [a]
uniteElems _ [] bs = bs
uniteElems _ as [] = as
uniteElems f (a:as) (b:bs) = f a b:uniteElems f as bs

markSelected :: (Int, Int) -> GElem -> [[GElem]] -> [[GElem]]
markSelected (i, j) (G v _ _) as = changeElem i (changeElem j (G v True True) (as!!i)) as

markVertical :: Int -> [[GElem]] -> [[GElem]]
markVertical j as = map (\a -> markHelper j (a!!j) a) as
    where
    markHelper :: Int -> GElem -> [GElem] -> [GElem]
    markHelper i (G a marked _) as = changeElem i (G a marked True) as

markHorizontal :: Int -> [[GElem]] -> [[GElem]]
markHorizontal i as = changeElem i a as
    where
    a = map (\(G a marked _) -> (G a marked True)) (as!!i)

markNeighbours :: (Int, Int) -> GElem -> [[GElem]] -> [[GElem]]
markNeighbours (i, j) g as = markHorizontal i (markVertical j (markSelected (i, j) g as))

divideSelectedRow :: Int -> GElem -> [[GElem]] -> [[GElem]]
divideSelectedRow i a as = changeElem i (map (/a) (as!!i)) as

decrOtherRows :: (Int,Int) -> GElem -> [[GElem]] -> [[GElem]]
decrOtherRows (i,j) g as = addElemAt i gs (decrHelper (removeElemAt i as)) 
    where
    gs = as!!i
    decrHelper :: [[GElem]] -> [[GElem]]
    decrHelper ls = map mapHelper ls
     where
     mapHelper l = uniteElems uniteHelper l gs
      where
      uniteHelper e1 e2 = e1 - ((l!!j) / g) * e2

hasSameLengthLists :: [[a]] -> Bool
hasSameLengthLists (a:as) = hasHelper a as
    where
    hasHelper :: [a] -> [[a]] -> Bool
    hasHelper _ [] = True
    hasHelper a (l:ls)
     | length a == length l = hasHelper a ls
     | otherwise = False

calcGauss :: [[Rational]] -> [[GElem]]
calcGauss [] = error "empty list"
calcGauss as
 | hasSameLengthLists as = gaussHelper1 $ map toGList as
 | otherwise = error "varying list lengths"

gaussHelper1 :: [[GElem]] -> [[GElem]]
gaussHelper1 as = gaussHelper2 (markVertical (length (as!!0) - 1) as)

gaussHelper2 :: [[GElem]] -> [[GElem]]
gaussHelper2 as
 | found = gaussHelper2 (decrOtherRows (i,j) (G 1 True True) (divideSelectedRow i g (markNeighbours (i,j) g (markSelected (i,j) g as))))
 | otherwise = as
    where
    (found, (i, j), g) = gaussSelector as

gauss :: [[Rational]] -> IO ()
gauss as = printGauss (calcGauss as)