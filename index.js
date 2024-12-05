/* Liquid Selector Variables List */
const toggleLiquidSelector = document.getElementById('liquidSelector');
const liquidForm = document.getElementById('liquidForm');
const liquidFieldset = liquidForm.querySelector('fieldset');
/* Gas Selector Variables List */
const toggleGasSelector = document.getElementById('gasSelector');
const gasForm = document.getElementById('gasForm');
const gasFieldset = gasForm.querySelector('fieldset');

/* Liquid Form Toggle Button */
toggleLiquidSelector.addEventListener('click', () => 
{
    const isToggled = toggleLiquidSelector.getAttribute('data-toggled') === 'true';
    toggleLiquidSelector.setAttribute('data-toggled', !isToggled); // Toggle the state

    /* Hide Gas Form */
    toggleGasSelector.setAttribute('data-toggled', 'false');
    gasForm.classList.add('hidden'); // Hide the form
    gasFieldset.disabled = !isToggled;

    // Update form and fieldset states
    liquidFieldset.disabled = isToggled;
    
    if(!isToggled)
    {
        liquidForm.classList.remove('hidden'); // Show the form
    }
    else
    {
        liquidForm.classList.add('hidden'); // Hide the form
    }
  
    //toggleLiquidSelector.textContent = isDisabled ? 'Enable Form' : 'Disable Form'; // Update button text
});

/* Gas Form Toggle Button */
toggleGasSelector.addEventListener('click', () => 
{
    const isToggled = toggleGasSelector.getAttribute('data-toggled') === 'true';
    toggleGasSelector.setAttribute('data-toggled', !isToggled); // Toggle the state
    
    /* Hide Liquid Form */
    toggleLiquidSelector.setAttribute('data-toggled', 'false');
    liquidForm.classList.add('hidden'); // Hide the form
    liquidFieldset.disabled = !isToggled;

    // Update form and fieldset states
    gasFieldset.disabled = isToggled;
    
    if(!isToggled)
    {
        gasForm.classList.remove('hidden'); // Show the form
    }
    else
    {
        gasForm.classList.add('hidden'); // Hide the form
    }
    
    //toggleGasSelector.textContent = isDisabled ? 'Enable Form' : 'Disable Form'; // Update button text
});


const liquidFluidSelect = document.getElementById('liquidFluidSelect');
//const actionButton = document.getElementById('calculate');

liquidFluidSelect.addEventListener('change', () => 
{
    // Update the button text with the selected option's text
    //actionButton.textContent = liquidFluidSelect.value;
});


/***********************************************************
const checkbox = document.getElementById('toggleCheckbox');
const content = document.getElementById('hiddenContent');
const extraField = document.getElementById('extraField');
const conditionalInput = document.getElementById('conditionalInput');

checkbox.addEventListener('change', () => 
{
    content.style.display = checkbox.checked ? 'block' : 'none';
});

checkbox.addEventListener('change', () => 
{
    if (checkbox.checked) 
    {

        extraField.classList.remove('hidden');
        conditionalInput.setAttribute('required', 'required');
    } 
    else 
    {
        extraField.classList.add('hidden');
        conditionalInput.removeAttribute('required');
    }
});
*/
