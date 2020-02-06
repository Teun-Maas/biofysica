classdef lsldert_cluster < lsldert_abstract_client
   properties
       clients    
   end
   
   methods
       function this = lsldert_cluster()
           this.clients = cell(0); 
       end
       
       function add_client(this, client)
           assert(isa(client,lsldert_client));
           n = numel(this.clients);
           this.clients{n+1} = client; 
       end
       
       function result=send(this, msg)
           n = numel(this.clients);
           result=cell(n,1);
           request=sprintf('%s\0',msg);
           
           % send request to all clients
           ii=1;
           while ii <= n
               this.clients{ii}.socket.send(request);
               ii=ii+1; 
           end

           % collect responses from all clients
           ii=1;
           while ii <= n
               rbytes=this.clients{ii}.socket.recv();
               result{ii}=reshape(char(rbytes),1,[]);
           end   
       end
   end    
end