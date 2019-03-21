// upload

slot_num = 1;
$('#upload_slot_add').on('click', function()
{
  slot_num += 1;
  $('#upload_slot').append('<li id="slot'+slot_num+'"><input type="file" name="file'+slot_num+'"></li>');
});
$('#upload_slot_remove').on('click', function()
{
  if (slot_num <= 1)
    return;

  $('#slot'+slot_num).remove();
  slot_num --;
});
