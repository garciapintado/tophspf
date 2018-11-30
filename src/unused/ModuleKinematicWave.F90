Module ModuleKinematicWave

  USE nrtype
  USE iso_varying_string
  USE strings, only : split


subroutine KinematicWave ()

        !Arguments-------------------------------------------------------------
        !real                                        :: LocalDT

        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        real                                        :: Slope
        real                                        :: level_left, level_right
        real                                        :: level_bottom, level_top
        real                                        :: HydraulicRadius, WetPerimeter
        real                                        :: Margin1, Margin2
        real                                        :: WaterDepth, MaxBottom
        integer                                     :: CHUNK, di, dj
        real(8)                                     :: MaxFlow
        !character(len=StringLength)                 :: Direction

        CHUNK = ChunkJ !CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !$OMP PARALLEL PRIVATE(I,J, Slope, level_left, level_right, level_bottom, level_top, &
        !$OMP HydraulicRadius, MaxFlow, Margin1, Margin2, WaterDepth, MaxBottom, WetPerimeter, di, dj)

        !X
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        xdo: do j = JLB, JUB
        ydo: do i = ILB, IUB
            if (Me%ComputeFaceU(i, j) == Compute) then

                !Adds to the final level the height of the buidings, if any
                if (Me%HydrodynamicApproximation == KinematicWave_) then

                    if (Me%Buildings) then
                        level_left  = Me%ExtVar%Topography(i, j-1) + Me%BuildingsHeight(i, j-1)
                        level_right = Me%ExtVar%Topography(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_left  = Me%ExtVar%Topography(i, j-1)
                        level_right = Me%ExtVar%Topography(i, j)
                    endif

                else if (Me%HydrodynamicApproximation == DiffusionWave_) then

                    if (Me%Buildings) then
                        level_left  = Me%myWaterLevel(i, j-1) + Me%BuildingsHeight(i, j-1)
                        level_right = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_left  = Me%myWaterLevel(i, j-1)
                        level_right = Me%myWaterLevel(i, j)
                    endif

                else

                    write(*,*)'Internal error'

                endif

                !Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_left - level_right) / Me%ExtVar%DZX(i, j-1))
                else
                    Slope           = (level_left - level_right) / Me%ExtVar%DZX(i, j-1)
                endif

                !Hydraulic Radius
!                Direction = "X"
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_left,level_right)
                !Wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DYY(i, j)

                !only compute in water column as MaxBottom (topography stairs descritization)
                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Water Depth consistent with AreaU computed (only water above max bottom)
                    WaterDepth = Me%AreaU(i,j) / Me%ExtVar%DYY(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i, j-1))

                    !to check wich cell to use to use since areaU depends on higher water level and max bottom
                    if (level_left .gt. level_right) then
                        dj = -1
                    else
                        dj = 0
                    endif

                    !Bottom Difference to adjacent cells (to check existence of margins on the side)
                    Margin1 = Me%ExtVar%Topography(i+1, j + dj) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i-1, j + dj) - MaxBottom

                    !if positive than there is a margin on the side and friction occurs at wet length
                    !If not basin points than result will be negative.
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif

                HydraulicRadius = Me%AreaU(i, j) / WetPerimeter


                !
                !MANNING'S EQUATION -  KINEMATIC WAVE
                !
                !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                if (Slope >= 0.0) then
                    Me%lFlowX(i, j) = Me%AreaU(i, j) * HydraulicRadius**(2./3.) * sqrt(Slope)          &
                                      / Me%OverlandCoefficientX(i,j)
                else
                    Me%lFlowX(i, j) = - Me%AreaU(i, j) * HydraulicRadius**(2./3.) * sqrt(-1.0 * Slope) &
                                      / Me%OverlandCoefficientX(i,j)
                endif


                !Limits Velocity to celerity if a free drop exists
                if (Me%HydrodynamicApproximation == DiffusionWave_ .and. Me%LimitToCriticalFlow) then
                    if ((level_left .lt. Me%ExtVar%Topography(i,j)) .or. (level_right .lt. Me%ExtVar%Topography(i,j-1))) then

                        !already defined in shorter
                        !WaterDepth = max (level_left, level_right) - max(Me%ExtVar%Topography(i, j-1), Me%ExtVar%Topography(i, j))
                        WaterDepth      = Me%AreaU(i, j)/Me%ExtVar%DYY(i,j)
                        MaxFlow         = Me%AreaU(i, j) * sqrt(Gravity * WaterDepth)
                        Me%lFlowX(i, j) = Min (MaxFlow, Me%lFlowX(i, j))

                    endif

                endif

            else

                Me%lFlowX(i, j) = 0.0

            endif

        enddo ydo
        enddo xdo
        !$OMP END DO NOWAIT

        !Y
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ComputeFaceV(i, j) == Compute) then

                !Adds to the final level the height of the buidings, if any
                if (Me%HydrodynamicApproximation == KinematicWave_) then

                    !Adds to the final level the height of the buidings, if any
                    if (Me%Buildings) then
                        level_bottom = Me%ExtVar%Topography(i-1, j) + Me%BuildingsHeight(i-1, j)
                        level_top    = Me%ExtVar%Topography(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_bottom = Me%ExtVar%Topography(i-1, j)
                        level_top    = Me%ExtVar%Topography(i, j)
                    endif

                else if (Me%HydrodynamicApproximation == DiffusionWave_) then

                    !Adds to the final level the height of the buidings, if any
                    if (Me%Buildings) then
                        level_bottom = Me%myWaterLevel(i-1, j) + Me%BuildingsHeight(i-1, j)
                        level_top    = Me%myWaterLevel(i, j)   + Me%BuildingsHeight(i, j  )
                    else
                        level_bottom = Me%myWaterLevel(i-1, j)
                        level_top    = Me%myWaterLevel(i, j)
                    endif

                else

                    write(*,*)'Internal error'

                endif


                !Slope
                if (Me%AdjustSlope) then
                    Slope           = AdjustSlope((level_bottom - level_top) / Me%ExtVar%DZY(i-1, j))
                else
                    Slope           = (level_bottom - level_top) / Me%ExtVar%DZY(i-1, j)
                endif

                !Hydraulic Radius
!                Direction = "Y"
!               !This function produced an overhead in openmp and the simulation took
!               !double the time so it was abandoned
!                HydraulicRadius = HydraulicRadius(i,j,Direction,level_bottom,level_top)
                !Wet perimeter, first is bottom
                WetPerimeter = Me%ExtVar%DXX(i, j)

                if ((Me%FaceWaterColumn == WCMaxBottom_) .and. (Me%CalculateCellMargins)) then
                    !Water Depth consistent with AreaV computed (only water above max bottom)
                    WaterDepth = Me%AreaV(i,j) / Me%ExtVar%DXX(i, j)
                    MaxBottom = max(Me%ExtVar%Topography(i, j), Me%ExtVar%Topography(i-1, j))

                    !to check wich cell to use since areaV depends on higher water level
                    if (level_bottom .gt. level_top) then
                        di = -1
                    else
                        di = 0
                    endif

                    !Bottom Difference to adjacent cells (to check existence of margins on the side)
                    Margin1 = Me%ExtVar%Topography(i + di,j+1) - MaxBottom
                    Margin2 = Me%ExtVar%Topography(i + di,j-1) - MaxBottom

                    !if positive than there is a margin on the side and friction occurs at wet length
                    if (Margin1 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin1)
                    endif
                    if (Margin2 .gt. 0.0) then
                        WetPerimeter = WetPerimeter + min(WaterDepth, Margin2)
                    endif
                endif

                !m = m2 / m
                HydraulicRadius = Me%AreaV(i, j) / WetPerimeter


                !
                !MANNING'S EQUATION -  KINEMATIC WAVE
                !
                !m3.s-1 = m2 * m(2/3) / (s.m(-1/3)) = m(8/3) * m(1/3) / s = m3.s-1
                if (Slope >= 0.0) then
                    Me%lFlowY(i, j) = Me%AreaV(i, j) * HydraulicRadius**(2./3.) * sqrt(Slope)            &
                                      / Me%OverlandCoefficientY(i,j)
                else
                    Me%lFlowY(i, j) = - Me%AreaV(i, j) * HydraulicRadius**(2./3.) * sqrt(-1.0 * Slope)   &
                                      / Me%OverlandCoefficientY(i,j)
                endif

                !Limits Velocity to reasonable values
                if (Me%HydrodynamicApproximation == DiffusionWave_ .and. Me%LimitToCriticalFlow) then

                    if ((level_bottom .lt. Me%ExtVar%Topography(i,j)) .or. (level_top .lt. Me%ExtVar%Topography(i-1,j))) then

                        !already defined in shorter
                        !WaterDepth = max (level_bottom, level_top) - max(Me%ExtVar%Topography(i-1,j), Me%ExtVar%Topography(i, j))
                        WaterDepth      = Me%AreaV(i, j)/Me%ExtVar%DXX(i,j)
                        MaxFlow         = Me%AreaV(i, j) * sqrt(Gravity * WaterDepth)
                        Me%lFlowY(i, j) = Min (MaxFlow, Me%lFlowY(i, j))

                    endif



                endif

            else

                Me%lFlowY(i, j) = 0.0

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine KinematicWave
