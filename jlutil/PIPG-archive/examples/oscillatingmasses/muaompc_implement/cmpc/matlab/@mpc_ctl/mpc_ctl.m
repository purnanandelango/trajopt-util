classdef mpc_ctl < handle
    properties(SetAccess=private)
        data
    end

    properties(Dependent)
        conf
        sys
        qpx
        wmx
        u_opt
        l_opt
        u_ref
        x_ref
    end
    
    methods(Static, Access=private)
        ctl = mpc_ctl_solve_problem(ctl, x)
        ctl = mpc_ctl_form_qp(ctl, x)
        ctl = mpc_ctl_get_data
    end

    methods
        function self = mpc_ctl()
            self.data = self.mpc_ctl_get_data;
        end

        function solve_problem(self, x)
            self.data = self.mpc_ctl_solve_problem(self.data, x);
        end
        
        function form_qp(self, x)
            self.data = self.mpc_ctl_form_qp(self.data, x);
        end

        function conf = get.conf(self)
            conf = self.data.conf;
        end

        function set.conf(self, conf)
            self.data.conf = conf;
        end

        function qpx = get.qpx(self)
            qpx = self.data.qpx;
        end

        function set.qpx(self, qpx)
            self.data.qpx = qpx;
        end

        function sys = get.sys(self)
            sys = self.data.sys;
        end

        function set.sys(self, sys)
            self.data.sys = sys;
        end

        function wmx = get.wmx(self)
            wmx = self.data.wmx;
        end

        function set.wmx(self, wmx)
            self.data.wmx = wmx;
        end

        function u_opt = get.u_opt(self)
            u_opt = self.data.u_opt;
        end

        function set.u_opt(self, u_opt)
            self.data.u_opt = u_opt;
        end

        function l_opt = get.l_opt(self)
            l_opt = self.data.l_opt;
        end

        function set.l_opt(self, l_opt)
            self.data.l_opt = l_opt;
        end

        function u_ref = get.u_ref(self)
            u_ref = self.data.u_ref;
        end

        function set.u_ref(self, u_ref)
            self.data.u_ref = u_ref;
        end

        function x_ref = get.x_ref(self)
            x_ref = self.data.x_ref;
        end

        function set.x_ref(self, x_ref)
            self.data.x_ref = x_ref;
        end
    end
end
