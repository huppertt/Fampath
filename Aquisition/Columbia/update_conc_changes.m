classdef update_conc_changes < handle
    properties
        HBO2=0;
        HHB=0;
        CCO=0;
        DPF;
        optode_dist;
        raw0;
        use_kalman=0;
        Q=1E-2;
    end
    properties(Hidden=true)
        lstValid_wavelengths;
        ext_coefs;
        inv_ext_coefs;
        
        modelKF;
        
        
    end
    
    methods
        function obj = update_conc_changes(wavelengths)
            [ext] = specific_extinction_coeffs_770to906;
            obj.ext_coefs=zeros(length(wavelengths),5);
            obj.ext_coefs(:,1)=wavelengths;
            for i=2:5
                obj.ext_coefs(:,i)=interp1(ext(:,1),ext(:,i),wavelengths);
            end
            obj.lstValid_wavelengths=find(~isnan(obj.ext_coefs(:,2)));
            obj.inv_ext_coefs = pinv(obj.ext_coefs(obj.lstValid_wavelengths,2:4)'*...
                                obj.ext_coefs(obj.lstValid_wavelengths,2:4))*...
                                obj.ext_coefs(obj.lstValid_wavelengths,2:4)';
            
            obj.modelKF=nirs.realtime.util.KalmanFilter(obj.Q*eye(3));             
                            
        end
        
        function set.Q(obj,Q)
            obj.Q=Q;
            if(isempty(obj.modelKF))
                obj.modelKF=nirs.realtime.util.KalmanFilter(obj.Q*eye(3));
            else
                obj.modelKF.Q=obj.Q*eye(3);
            end
        end
        
        function update(obj,raw)
            dOD = -log10(abs(raw(obj.lstValid_wavelengths)./obj.raw0(obj.lstValid_wavelengths)));
            dOD = dOD./obj.ext_coefs(obj.lstValid_wavelengths,5);
            
            if(obj.use_kalman)
                obj.modelKF.update(dOD,obj.ext_coefs(obj.lstValid_wavelengths,2:4)*(obj.optode_dist*obj.DPF))
                conc=obj.modelKF.B;
            else
                conc = obj.inv_ext_coefs*dOD/(obj.optode_dist*obj.DPF);
            end
            obj.HBO2=conc(1)*1000;
            obj.HHB =conc(2)*1000;
            obj.CCO =conc(3)*1000;
            
        end
        
    end
end