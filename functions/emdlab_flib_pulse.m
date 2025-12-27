function y = emdlab_flib_pulse(t,period,td,duty)

if duty>1, error('duty cannot be higher than 1.'); end
dwell = duty * period;

t = rem(t,period);
t1 = td;
t2 = td + dwell;
t2 = rem(t2,period);

y = 0;
if t1<=t2
    if (t>t1)&&(t<t2)
        y = 1;
    end
else
    if (t>t1)||(t<t2)
        y = 1;
    end
end

end