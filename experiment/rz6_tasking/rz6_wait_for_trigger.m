function desc=rz6_wait_for_trigger(trigger_type,external_trigger)

   desc = rz6_tq.define_task(rz6_tq.wait_for_trigger,0,0,trigger_type,external_trigger,0,0);

end
