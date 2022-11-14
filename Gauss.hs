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
    showList (xs) = (((intercalate " " $ map show (init xs)) ++ (" | " ++ show (last xs)))++)

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

findWhere :: (a -> Bool) -> [a] -> (Bool, Int, a)
findWhere _ [] = error "empty list"
findWhere f l@(a:as) = findHelper f a 0 l
    where
    findHelper ::(a -> Bool) ->  a -> Int -> [a] -> (Bool, Int, a)
    findHelper _ val _ [] = (False, 0, val)
    findHelper f val currInd (a:as)
     | f a = (True, currInd, a)
     | otherwise = findHelper f val (currInd+1) as

findWhereMulti :: (Ord a) => (a -> Bool) -> [[a]] -> (Bool, (Int, Int), a)
findWhereMulti f as = (found, newInd, val)
    where
    (found, ind, val) = findWhere f (joinLists as)
    newInd = (div ind (length (as !! 0)), mod ind (length (as !! 0)))

gaussSelector :: [[GElem]] -> (Bool, (Int, Int), GElem)
gaussSelector as = (found, (x, y), a)
 where (found, (x, y), a) = findWhereMulti (\(G a _ mNeigh) -> not mNeigh && a /= 0) as

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