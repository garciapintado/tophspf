! capabilities to include in tophspf



  ! simplified version of R sp SpatialPointsDataFrame class
  type T_SpatialPointsDataFrame
    TYPE (T_dataframe),      pointer                   :: data   => null()  ! attribute data as T_dataframe
    REAL,    DIMENSION(:,:), pointer                   :: coords => null()  ! coordinates matrix (points are rows in the matrix)
  end type T_SpatialPointsDataFrame

  TYPE T_Outlets





  INTEGER () outlets%code  :: vector with outlets
  outlets%hru
  outlets%

