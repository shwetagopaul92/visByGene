<script type='text/javascript'>
function ScrollDiv(){

   if(document.getElementById('ecran').scrollTop<(document.getElementById('ecran').scrollHeight-document.getElementById('ecran').offsetHeight)){-1
         document.getElementById('ecran').scrollTop=document.getElementById('ecran').scrollTop+1
         }
   else {document.getElementById('ecran').scrollTop=0;}
}

setInterval(ScrollDiv,50)
</script>
