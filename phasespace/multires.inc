
c variables for MultiResonances

      integer N_res(1:max_multires)
      double precision M2_X(1:max_multires), MG_X(1:max_multires), F_X(1:max_multires)
      double precision M2_i(1:max_res, 1:max_multires), MG_i(1:max_res,1:max_multires)
      double precision F_i(1:max_res, 1:max_multires), F_int(0:max_res, 1:max_multires)
      double precision lower(1:max_res, 1:max_multires), upper(1:max_res, 1:max_multires)
      common /PSmultires/ M2_X, MG_X, F_X, M2_i, MG_i, F_i, F_int, lower, upper, 
     &                    N_res

